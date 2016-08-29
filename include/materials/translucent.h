/****************************************************************************
 *    translucent.cc: translucent materials
 *    This is part of the yafray package
 *    Copyright (C) 2010  Ronnie
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

#ifndef TRANSLUCENT_H_
#define TRANSLUCENT_H_

// test review includes
#include <yafray_config.h>
#include <yafraycore/nodematerial.h>
#include <core_api/environment.h>
#include <core_api/scene.h>
#include <core_api/mcintegrator.h>
#include <materials/microfacet.h>



/*==========================
translucent material class
============================*/

__BEGIN_YAFRAY

#define C_TRANSLUCENT   0
#define C_GLOSSY        1
#define C_DIFFUSE       2

class translucentMat_t: public nodeMaterial_t
{
public:
    translucentMat_t(
        color_t diffuseC,   //! diffuse color
        color_t specC,      //! specular color
        color_t glossyC,    //! glossy color
        color_t siga,       //! sigma A color, reference to absorption or subsurface color
        color_t sigs,		//! sigma S color, reference to scatter color
        float sigs_factor,	//! sigma S factor
        float ior,			//! index of refraction
        float _g,			//! phase function
        float mT,			//! transmittance
        float mD,			//! diffuse reflect (diffusity)
        float mG,			//! glossy reflect (glossity)
        float exp);			//! fressnel exponent

    virtual ~translucentMat_t();

    virtual void initBSDF(const renderState_t &state, surfacePoint_t &sp, unsigned int &bsdfTypes)const;

    virtual color_t eval(const renderState_t &state, const surfacePoint_t &sp,
    		const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs)const;

    // povman: from material.h
    virtual color_t sample( const renderState_t &state, const surfacePoint_t &sp,
    		const vector3d_t &wo, vector3d_t &wi, sample_t &s, float &W)const;

    /*  UNUSED??
    virtual color_t sample(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo,
                      vector3d_t *const dir, color_t &tcol, sample_t &s, float *const W)const {return color_t(0.f);}
    */
    virtual color_t emit( const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const;

    virtual float pdf( const renderState_t &state, const surfacePoint_t &sp,
    		const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const;
    //virtual void getSpecular(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo,
    //                       bool &refl, bool &refr, vector3d_t *const dir, color_t *const col)const;
    static material_t* factory(paraMap_t &params, std::list< paraMap_t > &eparans, renderEnvironment_t &env);

protected:
    shaderNode_t* diffuseS;     //!< shader node for diffuse color
    shaderNode_t* glossyS;      //!< shader node for glossy color
    shaderNode_t* glossyRefS;   //!< shader node for glossy reflecity
    shaderNode_t* bumpS;        //!< shader node for bump mapping
    shaderNode_t *transpS;      //!< shader node for sigmaA (color_t)
    shaderNode_t *translS;      //!< shader node for sigmaS (color_t)

    color_t diffuseCol;
    color_t specRefCol;
    color_t gloss_color;

    float with_diffuse;
    float translucency, diffusity, glossity;
    float pDiffuse;
    float exponent;
    float sigmaS_Factor;

    BSDF_t cFlags[3];
    int nBSDF;

    // parameters for translucent property
    color_t sigma_s;
    color_t sigma_a;
    float   IOR;
    float   g;
};


__END_YAFRAY

#endif // TRANSLUCENT_H_
