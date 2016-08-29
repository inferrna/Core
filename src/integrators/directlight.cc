/****************************************************************************
*   directlight.cc: an integrator for direct lighting only
*   This is part of the yafaray package
*   Copyright (C) 2006  Mathias Wein (Lynx)
*   Copyright (C) 2009  Rodrigo Placencia (DarkTide)
*
*   This library is free software; you can redistribute it and/or
*   modify it under the terms of the GNU Lesser General Public
*   License as published by the Free Software Foundation; either
*   version 2.1 of the License, or (at your option) any later version.
*
*   This library is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*   Lesser General Public License for more details.
*
*   You should have received a copy of the GNU Lesser General Public
*   License along with this library; if not, write to the Free Software
*   Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/timer.h>
#include <sstream>
#include <iomanip>

__BEGIN_YAFRAY

class YAFRAYPLUGIN_EXPORT directLighting_t: public mcIntegrator_t
{
public:
    directLighting_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
    ~directLighting_t(); // povman: for SSS
    virtual bool preprocess();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
};

directLighting_t::directLighting_t(bool transpShad, int shadowDepth, int rayDepth)
{
    type = SURFACE;
    causRadius = 0.25;
    causDepth = 10;
    nCausPhotons = 100000;
    nCausSearch = 100;
    trShad = transpShad;
    usePhotonCaustics = false;
    sDepth = shadowDepth;
    rDepth = rayDepth;
    integratorName = "DirectLight";
    integratorShortName = "DL";
}

directLighting_t::~directLighting_t() // povman: for SSS
{
    destorySSSMaps();
}

bool directLighting_t::preprocess()
{
    bool success = true;
    std::stringstream set;
    settings = "";

    if(trShad)
    {
        set << "ShadowDepth: [" << sDepth << "]";
    }
    if(!set.str().empty()) set << "+";
    set << "RayDepth: [" << rDepth << "]";

    background = scene->getBackground();
    lights = scene->lights;

    if(usePhotonCaustics)
    {
        success = createCausticMap();
        if(!set.str().empty()) set << "+";
        set << "Caustics:" << nCausPhotons << " photons. ";
    }

    // povman: for SSS
    if (usePhotonSSS)
    {
        //success = createSSSMaps();
        success = createSSSMapsByPhotonTracing();
        set <<" SSS:["<< nSSSPhotons <<"] photons ";
        std::map<const object3d_t*, photonMap_t*>::iterator it = SSSMaps.begin();
        while (it!=SSSMaps.end())
        {
            it->second->updatePhTree();
            Y_INFO << integratorName << ": Builting: "  << it->second->nPhotons() << " photons. " << yendl;
            it++;
        }
    }
    // end

    if(useAmbientOcclusion)
    {
        if(!set.str().empty()) set << "+";
        set << "AO";
    }

    settings = set.str();

    return success;
}

colorA_t directLighting_t::integrate(renderState_t &state, diffRay_t &ray, colorPasses_t &colorPasses, int additionalDepth /*=0*/) const
{
    color_t col(0.0);
    CFLOAT alpha = 1.0;
    surfacePoint_t sp;
    void *o_udat = state.userdata;
    bool oldIncludeLights = state.includeLights;

	if (transpBackground) alpha = 0.0;

    // Shoot rays..
    if(scene->intersect(ray, sp))
    {
        unsigned char userdata[USER_DATA_SIZE];
        const material_t *material = sp.material;
        BSDF_t bsdfs;

        state.userdata = (void *) userdata;
        vector3d_t wo = -ray.dir;
        if(state.raylevel == 0) state.includeLights = true;

        material->initBSDF(state, sp, bsdfs);

        if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);

        if(bsdfs & BSDF_DIFFUSE)
        {
            col += estimateAllDirectLight(state, sp, wo);
            if (usePhotonCaustics)   col += estimateCausticPhotons(state, sp, wo);
            if (useAmbientOcclusion) col += sampleAmbientOcclusion(state, sp, wo);
        }

        // if have translucent SSS materials..
        if (bsdfs & BSDF_TRANSLUCENT)
        {
            // and use SSS photons is ON ..
            if (usePhotonSSS)
            {
                col += estimateSSSMaps(state, sp, wo);
                col += estimateSSSSingleSImportantSampling(state, sp, wo);
                //col += estimateSSSSingleScattering(state,sp,wo);
                //col += estimateSSSSingleScatteringPhotons(state,sp,wo);
            }
        }
        //
        recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);

        if(transpRefractedBackground)
        {
            CFLOAT m_alpha = material->getAlpha(state, sp, wo);
            alpha = m_alpha + (1.f-m_alpha)*alpha;
        }
        else alpha = 1.0;
    }
    else // Nothing hit, return background if any
    {
        if(background) col += (*background)(ray, state, false);
    }
	
    state.userdata = o_udat;
    state.includeLights = oldIncludeLights;

    //create new colorA_t for store also alpha value
    colorA_t result = col;
    colorA_t colVolTransmit = scene->volIntegrator->transmittance(state, ray);
    result *= colVolTransmit;
    result += scene->volIntegrator->integrate(state, ray);
	result.A = std::max(alpha, 1.f-colVolTransmit.A);
    
    return result;
}

integrator_t* directLighting_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
    bool transpShad=false;
    bool caustics=false;
    bool do_AO=false;
    int shadowDepth=5;
    int raydepth=5, cDepth=10;
    int search=100, photons=500000;
    int AO_samples = 32;
    double cRad = 0.25;
    double AO_dist = 1.0;
    color_t AO_col(1.f);
    bool bg_transp = true;
    bool bg_transp_refract = true;
    // SSS
    bool useSSS = false;
    int sssdepth = 10, sssPhotons = 100000;
    int singleSSamples = 128;
    float sScale = 40.f;

    params.getParam("raydepth", raydepth);
    params.getParam("transpShad", transpShad);
    params.getParam("shadowDepth", shadowDepth);
    params.getParam("caustics", caustics);
    params.getParam("photons", photons);
    params.getParam("caustic_mix", search);
    params.getParam("caustic_depth", cDepth);
    params.getParam("caustic_radius", cRad);
    params.getParam("do_AO", do_AO);
    params.getParam("AO_samples", AO_samples);
    params.getParam("AO_distance", AO_dist);
    params.getParam("AO_color", AO_col);
    params.getParam("bg_transp", bg_transp);
    params.getParam("bg_transp_refract", bg_transp_refract);
    //- SSS
    params.getParam("useSSS", useSSS);
    params.getParam("sssPhotons", sssPhotons);
    params.getParam("sssDepth", sssdepth);
    params.getParam("singleScatterSamples", singleSSamples);
    params.getParam("sssScale", sScale);
    //-

    directLighting_t *inte = new directLighting_t(transpShad, shadowDepth, raydepth);
    // caustic settings
    inte->usePhotonCaustics = caustics;
    inte->nCausPhotons = photons;
    inte->nCausSearch = search;
    inte->causDepth = cDepth;
    inte->causRadius = cRad;
    // AO settings
    inte->useAmbientOcclusion = do_AO;
    inte->aoSamples = AO_samples;
    inte->aoDist = AO_dist;
    inte->aoCol = AO_col;
    // Background settings
    inte->transpBackground = bg_transp;
    inte->transpRefractedBackground = bg_transp_refract;
    // sss settings
    inte->usePhotonSSS = useSSS;
    inte->nSSSPhotons = sssPhotons;
    inte->nSSSDepth = sssdepth;
    inte->nSingleScatterSamples = singleSSamples;
    inte->isDirectLight = true;
    inte->sssScale = sScale;

    //Y_INFO << "The translucent scale is " << inte->sssScale << yendl;
    // end
    return inte;
}

extern "C"
{

    YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
    {
        render.registerFactory("directlighting",directLighting_t::factory);
    }
}

__END_YAFRAY
