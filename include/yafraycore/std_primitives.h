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

#ifndef Y_STD_PRIMITIVE_H
#define Y_STD_PRIMITIVE_H

#include <core_api/primitive.h>

__BEGIN_YAFRAY

class renderEnvironment_t;
class paraMap_t;
class object3d_t;

class YAFRAYCORE_EXPORT sphere_t: public primitive_t
{
	public:
		sphere_t(point3d_t centr, PFLOAT rad, const material_t *m): center(centr), radius(rad), material(m) {}
		virtual bound_t getBound() const;
		virtual bool intersectsBound(exBound_t &b) const { return true; };
		//virtual bool clippingSupport() const { return false; }
		//virtual bool clipToBound(double bound[2][3], int axis, bound_t &clipped, void *d_old, void *d_new) const {return false;}
		virtual bool intersect(const ray_t &ray, PFLOAT *t, intersectData_t &data) const;
		virtual void getSurface(surfacePoint_t &sp, const point3d_t &hit, intersectData_t &data) const;
		virtual const material_t* getMaterial() const { return material; }
	protected:
		point3d_t center;
		PFLOAT radius;
		const material_t *material;
};

object3d_t* sphere_factory(paraMap_t &params, renderEnvironment_t &env);

__END_YAFRAY

#endif //Y_STD_PRIMITIVE_H
