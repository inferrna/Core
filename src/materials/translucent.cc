/****************************************************************************
 *    translucent.cc: translucent materials
 *    This is part of the yafray package
 *    Copyright (C) 2010  Ronnie
 *    Copyright (C) 2014  Pedro Alcaide
 *
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <materials/translucent.h>

__BEGIN_YAFRAY

translucentMat_t::translucentMat_t( color_t diffuseC, color_t specC, color_t glossyC, color_t siga, color_t sigs,
                                   float sigs_factor, float ior, float _g, float mT, float mD, float mG, float exp):
                                       diffuseS(0), glossyS(0), glossyRefS(0), bumpS(0), transpS(0), translS(0),
                                       diffuseCol(diffuseC), specRefCol(specC), gloss_color(glossyC), sigma_a(siga),
                                       sigma_s(sigs), sigmaS_Factor(sigs_factor), IOR(ior), g(_g), translucency(mT),
                                       diffusity(mD), glossity(mG), exponent(exp), pDiffuse(0.0)
{

    cFlags[C_TRANSLUCENT] = (BSDF_TRANSLUCENT);
    if (glossity > 0)
    {
        cFlags[C_GLOSSY] = (BSDF_GLOSSY | BSDF_REFLECT);

        if(diffusity > 0)
        {
            cFlags[C_DIFFUSE] = (BSDF_DIFFUSE | BSDF_REFLECT);
            with_diffuse = true;
            nBSDF = 3;
        }
        else
        {
            cFlags[C_DIFFUSE] = BSDF_NONE;
            nBSDF = 2;
        }
    }
    else
    {
        cFlags[C_GLOSSY] = cFlags[C_DIFFUSE] = BSDF_NONE;
        nBSDF = 1;
    }
    bsdfFlags = (cFlags[C_TRANSLUCENT] | cFlags[C_GLOSSY] | cFlags[C_DIFFUSE]);
}

translucentMat_t::~translucentMat_t()
{
    // empty
}
//-
void translucentMat_t::initBSDF(const renderState_t &state, surfacePoint_t &sp, unsigned int &bsdfTypes)const
{
    TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
    dat->stack = (char*)state.userdata + sizeof(TranslucentData_t);
    nodeStack_t stack(dat->stack);
    if(bumpS) evalBump(stack, state, sp, bumpS);

    std::vector<shaderNode_t *>::const_iterator iter, end=allViewindep.end();
    for(iter = allViewindep.begin(); iter!=end; ++iter) (*iter)->eval(stack, state, sp);

    dat->difC = diffuseS ? diffuseS->getColor(stack):diffuseCol;
    dat->sig_s = sigmaS_Factor * (translS ? translS->getColor(stack):this->sigma_s);
    dat->sig_a = transpS ? transpS->getColor(stack):this->sigma_a;
    dat->IOR = this->IOR;
    dat->g = this->g;

    dat->mDiffuse = this->diffusity;
    dat->mGlossy = glossyRefS ? glossyRefS->getScalar(stack) : this->glossity;
    dat->mTransl = this->translucency;

    dat->pDiffuse = std::min(0.6f , 1.f - (dat->mGlossy/(dat->mGlossy + (1.f - dat->mGlossy)*dat->mDiffuse)) );

    bsdfTypes = bsdfFlags;
}
//-
color_t translucentMat_t::eval(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs) const
{
    if( !(bsdfs & BSDF_DIFFUSE) || ((sp.Ng * wl) * (sp.Ng * wo)) < 0.f )
    {
        return color_t(0.f);
    }

    TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
    color_t col(0.f);
    bool diffuse_flag = bsdfs & BSDF_DIFFUSE;

    nodeStack_t stack(dat->stack);
    vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);

    float wiN = std::fabs(wl * N);
    float woN = std::fabs(wo * N);

    float Kr, Kt;
    fresnel(wl, N, IOR, Kr, Kt);

    float mR = (1.0f - Kt * dat->mTransl);

    // povman: if glossy_reflect value is > 0.0..
    if(bsdfs & BSDF_GLOSSY)
    {
        vector3d_t H = (wo + wl).normalize(); // half-angle
        float cos_wi_H = std::max(0.f, wl*H);
        float glossy;

        //glossy = Blinn_D(H*N, exponent) * SchlickFresnel(cos_wi_H, dat->mGlossy) / ASDivisor(cos_wi_H, woN, wiN);
        glossy = mR * Blinn_D(H * N, exponent) * SchlickFresnel(cos_wi_H, dat->mGlossy) / ASDivisor(cos_wi_H, woN, wiN);
        //-
        col = glossy * (glossyS ? glossyS->getColor(stack) : gloss_color); //col = glossy*gloss_color;
    }

    if(with_diffuse && diffuse_flag)
    {
        //col += diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diff_color)) * ((orenNayar)?OrenNayar(wi, wo, N):1.f);
        col += mR * diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diffuseCol));
    }

    return col;
}

color_t translucentMat_t::sample(const renderState_t &state, const surfacePoint_t &sp,
		const vector3d_t &wo, vector3d_t &wi, sample_t &samp, float &W)const
{
    TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
    float cos_Ng_wo = sp.Ng*wo;
    float cos_Ng_wi;
    vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
    vector3d_t Hs(0.f);
    samp.pdf = 0.f;
    float Kr, Kt;
    float wiN = 0.f , woN = 0.f;

    fresnel(wi, N, IOR, Kr, Kt);

    // missing! get components
    nodeStack_t stack(dat->stack);
    bool use[3] = {false, false, false};
    float sum = 0.f, accumC[3], val[3], width[3];
    int cIndex[3];	// entry values: 0 := specular part, 1 := glossy part, 2:= diffuse part;
    int rcIndex[3];	// reverse fmapping of cIndex, gives position of spec/glossy/diff in val/width array
    accumC[0] = Kt * dat->mTransl;
    accumC[1] = (1.f - Kt * dat->mTransl) * (1.f - dat->pDiffuse);
    accumC[2] = (1.f - Kt * dat->mTransl) * (dat->pDiffuse);

    int nMatch = 0, pick = -1;
    for(int i = 0; i < nBSDF; ++i)
    {
        if((samp.flags & cFlags[i]) == cFlags[i])
        {
            use[i] = true;
            width[nMatch] = accumC[i];
            cIndex[nMatch] = i;
            rcIndex[i] = nMatch;
            sum += width[nMatch];
            val[nMatch] = sum;
            ++nMatch;
        }
    }
    //-
    if(!nMatch || sum < 0.00001)
    {
        return color_t(0.f);
    }
    else if(nMatch==1)
    {
        pick=0;
        width[0]=1.f;
    }
    else
    {
        float inv_sum = 1.f/sum;
        for(int i=0; i<nMatch; ++i)
        {
            val[i] *= inv_sum;
            width[i] *= inv_sum;
            if((samp.s1 <= val[i]) && (pick < 0)){
                pick = i;
            }
        }
    }
    //-
    if(pick<0) pick=nMatch-1;
    //-
    float s1;
    if(pick>0)
    {
        s1 = (samp.s1 - val[pick-1]) / width[pick];
    }
    else
    {
        s1 = samp.s1 / width[pick];
    }

    color_t scolor(0.f);
	switch(cIndex[pick])
	{
    case C_TRANSLUCENT:
        // specular reflect
		break;
    case C_GLOSSY:
		// glossy; compute sampled half-angle vector H for Blinn distribution (microfacet.h)
		Blinn_Sample(Hs, s1, samp.s2, exponent);
		break;
	case C_DIFFUSE:
        // lambertian
		default:
		//! Sample a cosine-weighted hemisphere given the coordinate system built by N, Ru, Rv.(sample_utils.h)
		wi = SampleCosHemisphere(N, sp.NU, sp.NV, s1, samp.s2);
		cos_Ng_wi = sp.Ng*wi;
		if(cos_Ng_wo*cos_Ng_wi < 0) return color_t(0.f);
		break;
    }

    wiN = std::fabs(wi * N);
    woN = std::fabs(wo * N);

    if(cIndex[pick] != C_TRANSLUCENT)
    {
        //-------------------------------
        // if glossy reflection is > 0.0
        //-------------------------------
        if(use[C_GLOSSY])
        {
            float glossy;	//PFLOAT glossy;
            float cos_wo_H;//PFLOAT cos_wo_H;
            if(cIndex[pick] != C_GLOSSY)
            {
                vector3d_t H = (wi+wo).normalize();
                Hs = vector3d_t(H*sp.NU, H*sp.NV, H*N);
                cos_wo_H = wo*H;
            }
            else
            {
                vector3d_t H = Hs.x*sp.NU + Hs.y*sp.NV + Hs.z*N;
                cos_wo_H = wo*H;
                if ( cos_wo_H < 0.f )
                {
                    H.reflect(N);
                    cos_wo_H = wo*H;
                }
                // Compute incident direction by reflecting wo about H
                wi = reflect_dir(H, wo);
                cos_Ng_wi = sp.Ng*wi;
                if(cos_Ng_wo*cos_Ng_wi < 0) return color_t(0.f);
            }

            wiN = std::fabs(wi * N);

            samp.pdf += Blinn_Pdf(Hs.z, cos_wo_H, exponent) * width[rcIndex[C_GLOSSY]];
            glossy = Blinn_D(Hs.z, exponent) * SchlickFresnel(cos_wo_H, dat->mGlossy) / ASDivisor(cos_wo_H, woN, wiN);
            scolor = (float)glossy * (1.f-Kt * dat->mTransl) * (glossyS ? glossyS->getColor(stack) : gloss_color);
        }
        //---------------------------------
        // If 'diffuse reflection' is > 0.0
        //---------------------------------
        if(use[C_DIFFUSE])
        {
            scolor += (1.f-Kt*dat->mTransl)*diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diffuseCol));
            samp.pdf += wiN * width[rcIndex[C_DIFFUSE]];
        }
    }
    samp.sampledFlags = cFlags[cIndex[pick]];

    return scolor;
}

color_t translucentMat_t::emit(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const
{
    return color_t(0.f);
}

float translucentMat_t::pdf(const renderState_t &state, const surfacePoint_t &sp,
		const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const
{
    TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
    if( ((sp.Ng * wo) * (sp.Ng * wi)) < 0.f ) return 0.f;
    //-
    vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
    float pdf = 0.f;
    CFLOAT Kr, Kt;

    fresnel(wi, N, IOR, Kr, Kt);

    float accumC[3], sum=0.f, width;
    accumC[0] = Kt * dat->mTransl;
    accumC[1] = (1.f - Kt * dat->mTransl) * (1.f - dat->pDiffuse);
    accumC[2] = (1.f - Kt * dat->mTransl) * (dat->pDiffuse);

    int nMatch=0;
    for(int i=0; i<nBSDF; ++i)
    {
        if((bsdfs & cFlags[i]) == cFlags[i])
        {
            width = accumC[i];
            sum += width;
            if(i == C_GLOSSY)
            {
                vector3d_t H = (wi + wo).normalize();
                PFLOAT cos_wo_H = wo * H;
                PFLOAT cos_N_H = N * H;
                pdf += Blinn_Pdf(cos_N_H, cos_wo_H, exponent) * width;
            }
            else if(i == C_DIFFUSE)
            {
                pdf += std::fabs(wi * N) * width;
            }
            ++nMatch;
        }
    }
    if(!nMatch || sum < 0.00001) return 0.f;

    return pdf / sum;
}

//void translucentMat_t::getSpecular(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo,
//                         bool &refl, bool &refr, vector3d_t *const dir, color_t *const col)const
//{
//    PFLOAT cos_Ng_wo = sp.Ng*wo, cos_Ng_wi;
//
//    vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
//
//    float Kr, Kt;
//    fresnel(wo, N, IOR, Kr, Kt);
//
//    refr = false;
//    dir[0] = reflect_plane(N, wo);
//    //col[0] = (mirColS ? mirColS->getColor(stack) : specRefCol) * Kr;
//    col[0] = specRefCol * Kr;
//    refl = true;
//}

material_t* translucentMat_t::factory(paraMap_t &params, std::list< paraMap_t > &paramList, renderEnvironment_t &render)
{
    color_t col(1.0f);
    color_t glossyC(1.0f);
    color_t specC(1.0f);
    color_t sigA(0.01f);
    color_t sigS(1.0f);
    float sigS_factor = 1.0f;
    float ior = 1.3;
    float _g = 0;
    float mT = 0.9, mG = 1.0, mD = 0.001f;
    float expn = 800;

    params.getParam("color", col);
    params.getParam("glossy_color", glossyC);
    params.getParam("specular_color", specC);
    params.getParam("sigmaA", sigA);
    params.getParam("sigmaS", sigS);
    params.getParam("sigmaS_factor", sigS_factor);
    params.getParam("IOR", ior);
    params.getParam("g", _g);
    params.getParam("diffuse_reflect", mD);
    params.getParam("glossy_reflect", mG);
    params.getParam("sss_transmit", mT);
    params.getParam("exponent", expn);

    translucentMat_t *mat = new translucentMat_t(col, specC, glossyC, sigA, sigS, sigS_factor, ior, _g, mT, mD, mG, expn);

    std::vector<shaderNode_t *> roots;
    std::map<std::string, shaderNode_t *> nodeList;
    std::map<std::string, shaderNode_t *>::iterator actNode;

    // Prepare our node list
    nodeList["diffuse_shader"] = NULL;
    nodeList["glossy_shader"] = NULL;
    nodeList["glossy_reflect_shader"] = NULL;
    nodeList["bump_shader"] = NULL;
    nodeList["sigmaS_shader"] = NULL;
    nodeList["sigmaA_shader"] = NULL;

    // parse nodes
    if(mat->loadNodes(paramList, render))
    {
        mat->parseNodes(params, roots, nodeList);
    }
    else Y_ERROR << "Translucent: loadNodes() failed!" << yendl;

    mat->diffuseS = nodeList["diffuse_shader"];
    mat->glossyS = nodeList["glossy_shader"];
    mat->glossyRefS = nodeList["glossy_reflect_shader"];
    mat->bumpS = nodeList["bump_shader"];
    mat->translS = nodeList["sigmaS_shader"];
    mat->transpS = nodeList["sigmaA_shader"];

    // solve nodes order
    if(!roots.empty())
    {
        std::vector<shaderNode_t *> colorNodes;

        mat->solveNodesOrder(roots);

        if(mat->diffuseS)   mat->getNodeList(mat->diffuseS, colorNodes);
        if(mat->glossyS)    mat->getNodeList(mat->glossyS, colorNodes);
        if(mat->glossyRefS) mat->getNodeList(mat->glossyRefS, colorNodes);
        if(mat->transpS)    mat->getNodeList(mat->transpS, colorNodes);
        if(mat->translS)    mat->getNodeList(mat->translS, colorNodes);
        mat->filterNodes(colorNodes, mat->allViewdep, VIEW_DEP);
        mat->filterNodes(colorNodes, mat->allViewindep, VIEW_INDEP);
        if(mat->bumpS)      mat->getNodeList(mat->bumpS, mat->bumpNodes);
    }
    mat->reqMem = mat->reqNodeMem + sizeof(TranslucentData_t);

    return mat;
}

extern "C"
{
    YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
    {
        render.registerFactory("translucent", translucentMat_t::factory);
    }
}

__END_YAFRAY
