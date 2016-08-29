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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <yafray_constants.h>

#include <core_api/vector3d.h>
#include <core_api/ray.h>

__BEGIN_YAFRAY

struct Plane
{
    vector3d_t p;
    vector3d_t n;
};

inline float ray_plane_intersection(ray_t const& ray, Plane const& plane)
{
    return plane.n * (plane.p - vector3d_t(ray.from)) / (ray.dir * plane.n);
}

__END_YAFRAY

#endif // GEOMETRY_H
