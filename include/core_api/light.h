/****************************************************************************
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef Y_LIGHT_H
#define Y_LIGHT_H

#include "ray.h"
#include "scene.h"

__BEGIN_YAFRAY

struct surfacePoint_t;
class background_t;

enum { LIGHT_NONE = 0, LIGHT_DIRACDIR = 1, LIGHT_SINGULAR = 1<<1 }; // "LIGHT_DIRACDIR" *must* be same as "BSDF_SPECULAR" (material.h)!
typedef unsigned int LIGHTF_t;

struct lSample_t
{
	lSample_t(surfacePoint_t *s_p=0): sp(s_p) {}
	//<! 2d sample value for choosing a surface point on the light.
	float s1, s2;
	//<! 2d sample value for choosing an outgoing direction on the light (emitSample)
	float s3, s4; 
	//<! "standard" directional pdf from illuminated surface point for MC integration of direct lighting (illumSample)
	float pdf; 
	//<! probability density for generating this sample direction (emitSample)
	float dirPdf;
	//<! probability density for generating this sample point on light surface (emitSample)
	float areaPdf;
	//<! color of the generated sample
	color_t col;
	//<! flags of the sampled light source
	LIGHTF_t flags;
	//!< surface point on the light source, may only be complete enough to call other light methods with it!
	surfacePoint_t *sp;
};

class light_t
{
	public:
		//! allow for preprocessing when scene loading has finished
		virtual void init(scene_t &scene) {}

		//! total energy emitted during whole frame
		virtual color_t totalEnergy() const = 0;

		//! emit a photon
		virtual color_t emitPhoton(float s1, float s2, float s3, float s4, ray_t &ray, float &ipdf) const = 0;

		//! create a sample of light emission, similar to emitPhoton, just more suited for bidirectional methods
		/*! fill in s.dirPdf, s.areaPdf, s.col and s.flags, and s.sp if not nullptr */
		virtual color_t emitSample(vector3d_t &wo, lSample_t &s) const{return color_t(0.f);};

		//! indicate whether the light has a dirac delta distribution or not
		virtual bool diracLight() const = 0;

		//! illuminate a given surface point, generating sample s, fill in s.sp if not NULL; Set ray to test visibility by integrator
		/*! fill in s.pdf, s.col and s.flags */
		virtual bool illumSample(const surfacePoint_t &sp, lSample_t &s, ray_t &wi) const = 0;

		//! illuminate a given surface point; Set ray to test visibility by integrator. Only for dirac lights.
		/*!	return false only if no light is emitted towards sp, e.g. outside cone angle of spot light	*/
		virtual bool illuminate(const surfacePoint_t &sp, color_t &col, ray_t &wi) const = 0;

		//! indicate whether the light can intersect with a ray (by the intersect function)
		virtual bool canIntersect() const { return false; }

		//! intersect the light source with a ray, giving back distance, energy and 1/PDF
		virtual bool intersect(const ray_t &ray, float &t, color_t &col, float &ipdf) const { return false; }

		//! get the pdf for sampling the incoming direction wi at surface point sp (illumSample!)
		/*! this method requires an intersection point with the light (sp_light). Otherwise, use intersect() */
		virtual float illumPdf(const surfacePoint_t &sp, const surfacePoint_t &sp_light) const { return 0.f; }

		//! get the pdf values for sampling point sp on the light and outgoing direction wo when emitting energy (emitSample, NOT illumSample)
		/*! sp should've been generated from illumSample or emitSample, and may only be complete enough to call light functions! */
		virtual void emitPdf(const surfacePoint_t &sp, const vector3d_t &wo, float &areaPdf, float &dirPdf, float &cos_wo) const { areaPdf=0.f; dirPdf=0.f; }
		
		//! checks if the light can shoot caustic photons (photon map integrator)
		virtual bool shootsCausticP() const { return true;}
		
		//! checks if the light can shoot diffuse photons (photon map integrator)
		virtual bool shootsDiffuseP() const { return true;}
		
		//! (preferred) number of samples for direct lighting
		virtual int nSamples() const { return 8; }
		virtual ~light_t() {}
		
		//! This method must be called right after the factory is called on a background light or the light will fail
		virtual void setBackground(background_t *bg) { background = bg; }
		//! Enable/disable entire light source
		bool lightEnabled() const { return lLightEnabled;}
		bool castShadows() const { return lCastShadows; }
		//! checks if the light can shoot caustic photons (photonmap integrator)
		bool shootsCausticP() const { return lShootCaustic; }
		//! checks if the light can shoot diffuse photons (photonmap integrator)
		bool shootsDiffuseP() const { return lShootDiffuse; }
		//! checks if the light is a photon-only light (only shoots photons, not illuminating)
		bool photonOnly() const { return lPhotonOnly; }
		//! sets clampIntersect value to reduce noise at the expense of realism and inexact overall lighting
		void setClampIntersect(float clamp) { lClampIntersect = clamp; }

		light_t(): flags(LIGHT_NONE),lLightEnabled(true),lCastShadows(true),lShootCaustic(true),lShootDiffuse(true),lPhotonOnly(false) {}
		light_t(LIGHTF_t _flags): flags(_flags) {}
		LIGHTF_t getFlags() const { return flags; }

	protected:
		LIGHTF_t flags;
		background_t* background;
	    bool lLightEnabled; //!< enable/disable light
		bool lCastShadows; //!< enable/disable if the light should cast direct shadows
		bool lShootCaustic; //!<enable/disable if the light can shoot caustic photons (photonmap integrator)
		bool lShootDiffuse; //!<enable/disable if the light can shoot diffuse photons (photonmap integrator)
		bool lPhotonOnly; //!<enable/disable if the light is a photon-only light (only shoots photons, not illuminating)
		float lClampIntersect = 0.f;	//!<trick to reduce light sampling noise at the expense of realism and inexact overall light. 0.f disables clamping

};

__END_YAFRAY

#endif // Y_LIGHT_H
