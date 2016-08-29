/****************************************************************************
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
#ifndef BDPT_H
#define BDPT_H

#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/tiledintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/camera.h>
#include <core_api/imagefilm.h>
#include <integrators/integr_utils.h>
#include <utilities/mcqmc.h>
#include <sstream>
#include <algorithm>


__BEGIN_YAFRAY

/*  conventions:
    y_0 := point on light source
    z_0 := point on camera lens
    x_0...x_k := path vertices, while x_0...x_s-1 are y_0...y_s-1 and x_k...x_s are z_0...z_t-1
    so x_i <=> z_k-i, for i>=s
    */


#define MAX_PATH_LENGTH 32
#define MIN_PATH_LENGTH 3

#define _BIDIR_DEBUG 0
//#define _DO_LIGHTIMAGE 1

/*! class that holds some vertex y_i/z_i (depending on wether it is a light or camera path)
*/

class pathVertex_t
{
public:
    surfacePoint_t sp;  //!< surface point at which the path vertex lies
    BSDF_t flags;       //!< flags of the sampled BSDF component (not all components of the sp!)
    color_t alpha;      //!< cumulative subpath weight; note that y_i/z_i stores alpha_i+1 !
    color_t f_s;        //!< f(x_i-1, x_i, x_i+1), i.e. throughput from last to next path vertex
    vector3d_t wi, wo;  //!< sampled direction for next vertex (if available)
    float ds;           //!< squared distance between x_i-1 and x_i
    float G;            //!< geometric factor G(x_i-1, x_i), required for MIS
    float qi_wo;        //!< russian roulette probability for terminating path
    float qi_wi;        //!< russian roulette probability for terminating path when generating path in opposite direction
    float cos_wi, cos_wo;   //!< (absolute) cosine of the incoming (wi) and sampled (wo) path direction
    float pdf_wi, pdf_wo;   //!< the pdf for sampling wi from wo and wo from wi respectively
    void *userdata;     //!< user data of the material at sp (required for sampling and evaluating)
};

/*! vertices of a connected path going forward from light to eye;
    conventions: path vertices are named x_0...x_k, with k=s+t-1 again.
    x_0 lies on the light source, x_k on the camera */

struct pathEvalVert_t
{
    bool specular;      //!< indicate that the ingoing direction determines the outgoing one (and vice versa)
    union
    {
        float pdf_f;    //!< pdf of sampling forward direction (x_i->x_i+1) given the backward direction
        float pdf_A_k;  //!< in case of lense vertex we have area pdf here, there is no forward path segment
    };
    union
    {
        float pdf_b;    //!< pdf of sampling backward direction (x_i-1->x_i) given forward direction
        float pdf_A_0;  //!< in case of light vertex we have area pdf here, there is no backward path segment
    };
    float G;            //!< geometric term G(x_i-1, x_i)
};

class pathData_t
{
public:
    std::vector<pathVertex_t> lightPath, eyePath;
    std::vector<pathEvalVert_t> path;
    //pathCon_t pc;
    // additional information for current path connection:
    vector3d_t w_l_e;   //!< direction of edge from light to eye vertex, i.e. y_s to z_t
    color_t f_y, f_z;   //!< f for light and eye vertex that are connected
    PFLOAT u, v;        //!< current position on image plane
    float d_yz;         //!< distance between y_s to z_t
    const light_t *light; //!< the light source to which the current path is connected
    //float pdf_Ad_0; //!< pdf for direct lighting strategy
    float pdf_emit, pdf_illum;  //!< light pdfs required to calculate p1 for direct lighting strategy
    bool singularL; //!< true if light has zero area (point lights for example)
    int nPaths;     //!< number of paths that have been sampled (for current thread and image)
};

class YAFRAYPLUGIN_EXPORT biDirIntegrator_t : public tiledIntegrator_t
{
public:
    biDirIntegrator_t(bool transpShad=false, int shadowDepth=4);
    virtual ~biDirIntegrator_t(){};
    virtual bool preprocess();
    virtual void cleanup();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
protected:
    int createPath(renderState_t &state, ray_t &start, std::vector<pathVertex_t> &path, int maxLen) const;
    color_t evalPath(renderState_t &state, int s, int t, pathData_t &pd) const;
    color_t evalLPath(renderState_t &state, int t, pathData_t &pd, ray_t &lRay, const color_t &lcol) const;
    color_t evalPathE(renderState_t &state, int s, pathData_t &pd) const;
    bool connectPaths(renderState_t &state, int s, int t, pathData_t &pd) const;
    bool connectLPath(renderState_t &state, int t, pathData_t &pd, ray_t &lRay, color_t &lcol) const;
    bool connectPathE(renderState_t &state, int s, pathData_t &pd) const;
    //color_t estimateOneDirect(renderState_t &state, const surfacePoint_t &sp, vector3d_t wo, pathCon_t &pc)const;
    CFLOAT pathWeight(renderState_t &state, int s, int t, pathData_t &pd) const;
    CFLOAT pathWeight_0t(renderState_t &state, int t, pathData_t &pd) const;

    background_t *background;
    const camera_t *cam;
    bool trShad;
    bool use_bg;    //!< configuration; include background for GI
    bool ibl;       //!< configuration; use background light, if available
    bool include_bg;    //!< determined on precrocess;
    int sDepth, rDepth, bounces;

    // povman test ---------------
    bool do_lightImage;         //
    //---------------------------- end

    std::vector<light_t*> lights;
    //mutable std::vector<pathVertex_t> lightPath, eyePath;
    //mutable int nPaths;
    //mutable pathData_t pathData;
    mutable std::vector<pathData_t> threadData;
    pdf1D_t *lightPowerD;
    float fNumLights;
    std::map <const light_t*, CFLOAT> invLightPowerD;
    imageFilm_t *lightImage;
    // test
    bool transpBackground; //! Render background as transparent
    bool transpRefractedBackground; //! Render refractions of background as transparent
};


__END_YAFRAY

#endif
