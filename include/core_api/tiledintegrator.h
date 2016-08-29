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

#ifndef Y_TILEDINTEGRATOR_H
#define Y_TILEDINTEGRATOR_H

#include <core_api/integrator.h>
#include <core_api/imagesplitter.h>
#include <core_api/material.h>

__BEGIN_YAFRAY

class threadControl_t;

class YAFRAYCORE_EXPORT tiledIntegrator_t: public surfaceIntegrator_t
{
	public:
		/*! Rendering prepasses to precalc suff in case needed */
		virtual void preRender(); //!< Called before the render starts and after the minDepth and maxDepth are calculated
		virtual void prePass(int samples, int offset, bool adaptive); //!< Called before the proper rendering of all the tiles starts
		virtual void preTile(renderArea_t &a, int n_samples, int offset, bool adaptive, int threadID); //!< Called brfore each tile is rendered
		
		/*! do whatever is required to render the image; default implementation renders image in passes
		dividing each pass into tiles for multithreading. */
		virtual bool render(int numView, imageFilm_t *imageFilm);
		/*! render a pass; only required by the default implementation of render() */
		virtual bool renderPass(int numView, int samples, int offset, bool adaptive, int AA_pass_number);
		/*! render a tile; only required by default implementation of render() */
		virtual bool renderTile(int numView, renderArea_t &a, int n_samples, int offset, bool adaptive, int threadID, int AA_pass_number = 0);
		virtual void renderWorker(int mNumView, tiledIntegrator_t *integrator, scene_t *scene, imageFilm_t *imageFilm, threadControl_t *control, int threadID, int samples, int offset=0, bool adaptive=false, int AA_pass=0);
		
//		virtual void recursiveRaytrace(renderState_t &state, diffRay_t &ray, int rDepth, BSDF_t bsdfs, surfacePoint_t &sp, vector3d_t &wo, color_t &col, float &alpha) const;
		virtual void precalcDepths();
		virtual void generateCommonRenderPasses(colorPasses_t &colorPasses, renderState_t &state, const surfacePoint_t &sp) const; //!< Generates render passes common to all integrators
	
	protected:
		int AA_samples, AA_passes, AA_inc_samples;
		float iAA_passes; //!< Inverse of AA_passes used for depth map
		float AA_threshold;
		float AA_resampled_floor; //!< minimum amount of resampled pixels (% of the total pixels) below which we will automatically decrease the AA_threshold value for the next pass
		float AA_sample_multiplier_factor;
		float AA_light_sample_multiplier_factor;
		float AA_indirect_sample_multiplier_factor;
		bool AA_detect_color_noise;
		int AA_dark_detection_type;
		float AA_dark_threshold_factor;
		int AA_variance_edge_size;
		int AA_variance_pixels;
		float AA_clamp_samples;
		float AA_clamp_indirect;
		float AA_sample_multiplier;
		float AA_light_sample_multiplier;
		float AA_indirect_sample_multiplier;
		imageFilm_t *imageFilm;
		float maxDepth; //!< Inverse of max depth from camera within the scene boundaries
		float minDepth; //!< Distance between camera and the closest object on the scene
		bool diffRaysEnabled;	//!< Differential rays enabled/disabled - for future motion blur / interference features

};

struct threadControl_t
{
	threadControl_t() : finishedThreads(0) {}
	std::mutex m;
	std::condition_variable c; //!< condition variable to signal main thread
	std::vector<renderArea_t> areas; //!< area to be output to e.g. blender, if any
	volatile int finishedThreads; //!< number of finished threads, lock countCV when increasing/reading!
};

__END_YAFRAY

#endif // Y_TILEDINTEGRATOR_H
