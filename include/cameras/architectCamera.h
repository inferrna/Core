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

#ifndef Y_ARCHITECTCAMERA_H
#define Y_ARCHITECTCAMERA_H

#include <cameras/perspectiveCamera.h>

__BEGIN_YAFRAY

class paraMap_t;
class renderEnvironment_t;

class architectCam_t: public perspectiveCam_t
{
	public:
        architectCam_t(const point3d_t &pos, const point3d_t &look, const point3d_t &up,
                       int _resx, int _resy, PFLOAT aspect=1,
                       PFLOAT df=1, PFLOAT ap=0, PFLOAT dofd=0, bokehType bt=BK_DISK1, bkhBiasType bbt=BB_NONE, PFLOAT bro=0,
                       float const near_clip_distance = 0.0f, float const far_clip_distance = 1e6f);
		virtual ~architectCam_t();
		virtual void setAxis(const vector3d_t &vx, const vector3d_t &vy, const vector3d_t &vz);
		virtual point3d_t screenproject(const point3d_t &p) const;

		static camera_t* factory(paraMap_t &params, renderEnvironment_t &render);
};

__END_YAFRAY

#endif // Y_ARCHITECTCAMERA_H
