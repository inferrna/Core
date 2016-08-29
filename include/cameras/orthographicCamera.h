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

#ifndef Y_ORTHOGRAPHICCAMERA_H
#define Y_ORTHOGRAPHICCAMERA_H

#include <yafray_config.h>

#include <core_api/camera.h>

__BEGIN_YAFRAY

class paraMap_t;
class renderEnvironment_t;

class orthoCam_t: public camera_t
{
	public:
		orthoCam_t(const point3d_t &pos, const point3d_t &look, const point3d_t &up,
                   int _resx, int _resy, PFLOAT aspect, PFLOAT scale,
                   float const near_clip_distance = 0.0f, float const far_clip_distance = 1e6f);
		virtual void setAxis(const vector3d_t &vx, const vector3d_t &vy, const vector3d_t &vz);
		virtual ray_t shootRay(PFLOAT px, PFLOAT py, float lu, float lv, PFLOAT &wt) const;
		virtual point3d_t screenproject(const point3d_t &p) const;
		
		static camera_t* factory(paraMap_t &params, renderEnvironment_t &render);
	protected:
		PFLOAT scale;
		point3d_t pos;
};

__END_YAFRAY

#endif // Y_ORTHOGRAPHICCAMERA_H
