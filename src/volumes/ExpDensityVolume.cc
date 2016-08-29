/*************************************************************************
 *	This file is part of the yafray package
 *  Copyright (C) 2002  Alejandro Conty Estévez
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include <core_api/ray.h>
#include <core_api/color.h>
#include <core_api/volume.h>
#include <core_api/bound.h>
#include <core_api/volume.h>
#include <core_api/surface.h>
#include <core_api/environment.h>

__BEGIN_YAFRAY

struct renderState_t;
struct pSample_t;

class ExpDensityVolume : public DensityVolume
{
public:

	ExpDensityVolume(color_t sa, color_t ss, color_t le, float gg, point3d_t pmin, point3d_t pmax, int attgridScale, float aa, float bb) :
		DensityVolume(sa, ss, le, gg, pmin, pmax, attgridScale)
	{
		a = aa;
		b = bb;
		Y_INFO << "ExpDensityVolume vol: " << s_a << " " << s_s << " " << l_e << " " << a << " " << b << yendl;
	}

	virtual float Density(point3d_t p);

	static VolumeRegion* factory(paraMap_t &params, renderEnvironment_t &render);

protected:
	float a, b;
};

float ExpDensityVolume::Density(point3d_t p)
{
	float height = p.z - bBox.a.z;
	return a * fExp(-b * height);
}

VolumeRegion* ExpDensityVolume::factory(paraMap_t &params,renderEnvironment_t &render)
{
	float ss = .1f;
	float sa = .1f;
	float le = .0f;
	float g = .0f;
	float a = 1.f;
	float b = 1.f;
	float min[] = {0, 0, 0};
	float max[] = {0, 0, 0};
	int attSc = 1;

	params.getParam("sigma_s", ss);
	params.getParam("sigma_a", sa);
	params.getParam("l_e", le);
	params.getParam("g", g);
	params.getParam("a", a);
	params.getParam("b", b);
	params.getParam("minX", min[0]);
	params.getParam("minY", min[1]);
	params.getParam("minZ", min[2]);
	params.getParam("maxX", max[0]);
	params.getParam("maxY", max[1]);
	params.getParam("maxZ", max[2]);
	params.getParam("attgridScale", attSc);

	ExpDensityVolume *vol = new ExpDensityVolume(color_t(sa), color_t(ss), color_t(le), g,
						point3d_t(min[0], min[1], min[2]), point3d_t(max[0], max[1], max[2]), attSc, a, b);
	return vol;
}

extern "C"
{
	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("ExpDensityVolume", ExpDensityVolume::factory);
	}
}

__END_YAFRAY
