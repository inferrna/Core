/*************************************************************************
 *	This is part of TheBounty package
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
 */

 #include <yafray_config.h>

#include <core_api/ray.h>
#include <core_api/color.h>
#include <core_api/volume.h>
#include <core_api/bound.h>
#include <core_api/surface.h>
#include <core_api/texture.h>
#include <core_api/environment.h>
#include <utilities/mcqmc.h>

#include <fstream>
#include <cstdlib>

__BEGIN_YAFRAY

class renderState_t;
class pSample_t;

class GridVolume : public DensityVolume
{
public:

	GridVolume(color_t sa, color_t ss, color_t le, float gg, point3d_t pmin, point3d_t pmax, int attSc, std::string gridfile):
		DensityVolume(sa, ss, le, gg, pmin, pmax, attSc)
	{
		bBox = bound_t(pmin, pmax);
		s_a = sa;
		s_s = ss;
		l_e = le;
		g = gg;
		haveS_a = (s_a.energy() > 1e-4f);
		haveS_s = (s_s.energy() > 1e-4f);
		haveL_e = (l_e.energy() > 1e-4f);

		std::ifstream inputStream;

		inputStream.open(gridfile.c_str());
		if(!inputStream)
		{
			Y_ERROR << "GridVolume: Error opening density file" << yendl;
		}
		else
		{
			Y_INFO << "Processing density file: " << gridfile.c_str() << yendl;
		}

		inputStream.seekg(0, std::ios_base::beg);
		std::ifstream::pos_type begin_pos = inputStream.tellg();
		inputStream.seekg(0, std::ios_base::end);
		int fileSize = static_cast<int>(inputStream.tellg() - begin_pos);
		fileSize -= 6;
		inputStream.seekg(0, std::ios_base::beg);

		int dim[3];
		for (int i = 0; i < 3; ++i)
		{
			short i0 = 0, i1 = 0;
			inputStream.read( (char*)&i0, 1 );
			inputStream.read( (char*)&i1, 1 );
			Y_INFO << "GridVolume: " << i0 << " " << i1 << yendl;
			dim[i] = (((unsigned short)i0 << 8) | (unsigned short)i1);
		}
		int sizePerVoxel = fileSize / (dim[0] * dim[1] * dim[2]);

		Y_INFO << "GridVolume: "<< dim[0] <<" "<< dim[1] <<" "<< dim[2] <<" File size: "<< fileSize <<" "<< sizePerVoxel << yendl;

		sizeX = dim[0];
		sizeY = dim[1];
		sizeZ = dim[2];

		/* povman: change to sizeof(float*) from sizeof(float) for fix crash on x64 system's.
		 * atm work fine on Ubuntu 14.04 x64, but is need make more test with others OS
		 */
		grid = (float***)malloc(sizeX * sizeof(float*));

		for (int x = 0; x < sizeX; ++x)
		{
			grid[x] = (float**)malloc(sizeY * sizeof(float*));

			for (int y = 0; y < sizeY; ++y)
			{
				grid[x][y] = (float*)malloc(sizeZ * sizeof(float*));
			}
		}
		//--
		for (int z = 0; z < sizeZ; ++z)
		{
			for (int y = 0; y < sizeY; ++y)
			{
				for (int x = 0; x < sizeX; ++x)
				{
					int voxel = 0;
					inputStream.read( (char*)&voxel, 1 );
					grid[x][y][z] = voxel / 255.f;
					/* povman: code commented by author..
					float r = sizeX / 2.f;
					float r2 = r*r;
					float dist = sqrt((x-r)*(x-r) + (y-r)*(y-r) + (z-r)*(z-r));
					grid[x][y][z] = (dist < r) ? 1.f-dist/r : 0.0f;
					*/
				}
			}
		}
		Y_INFO << "GridVolume: Vol.[" << s_a << ", " << s_s << ", " << l_e << "]" << yendl;
	}

	~GridVolume();

	virtual float Density(point3d_t p);

	static VolumeRegion* factory(paraMap_t &params, renderEnvironment_t &render);

protected:
	float*** grid;
	int sizeX, sizeY, sizeZ;
};

// test for expand class
GridVolume::~GridVolume()
{
	Y_INFO << "GridVolume: Freeing grid data" << yendl;

	for (int x = 0; x < sizeX; ++x)
	{
		for (int y = 0; y < sizeY; ++y)
		{
			free(grid[x][y]);
		}
		free(grid[x]);
	}
	free(grid);
}

inline float min(float a, float b){	return (a > b) ? b : a;}

inline float max(float a, float b){	return (a < b) ? b : a;}


float GridVolume::Density(const point3d_t p)
{
	float x = (p.x - bBox.a.x) / bBox.longX() * sizeX - .5f;
	float y = (p.y - bBox.a.y) / bBox.longY() * sizeY - .5f;
	float z = (p.z - bBox.a.z) / bBox.longZ() * sizeZ - .5f;

	int x0 = max(0, floor(x));
	int y0 = max(0, floor(y));
	int z0 = max(0, floor(z));

	int x1 = min(sizeX - 1, ceil(x));
	int y1 = min(sizeY - 1, ceil(y));
	int z1 = min(sizeZ - 1, ceil(z));

	float xd = x - x0;
	float yd = y - y0;
	float zd = z - z0;

	float i1 = grid[x0][y0][z0] * (1-zd) + grid[x0][y0][z1] * zd;
	float i2 = grid[x0][y1][z0] * (1-zd) + grid[x0][y1][z1] * zd;
	float j1 = grid[x1][y0][z0] * (1-zd) + grid[x1][y0][z1] * zd;
	float j2 = grid[x1][y1][z0] * (1-zd) + grid[x1][y1][z1] * zd;

	float w1 = i1 * (1 - yd) + i2 * yd;
	float w2 = j1 * (1 - yd) + j2 * yd;

	float dens = w1 * (1 - xd) + w2 * xd;

	return dens;
}

VolumeRegion* GridVolume::factory(paraMap_t &params, renderEnvironment_t &render)
{
	float ss = .1f;
	float sa = .1f;
	float le = .0f;
	float g = .0f;
	float min[] = {0, 0, 0};
	float max[] = {0, 0, 0};
	int attSc = 1;
	std::string densityFile; // = NULL;

	params.getParam("sigma_s", ss);
	params.getParam("sigma_a", sa);
	params.getParam("l_e", le);
	params.getParam("g", g);
	params.getParam("minX", min[0]);
	params.getParam("minY", min[1]);
	params.getParam("minZ", min[2]);
	params.getParam("maxX", max[0]);
	params.getParam("maxY", max[1]);
	params.getParam("maxZ", max[2]);
	params.getParam("attgridScale", attSc);
	params.getParam("density_file", densityFile);

	GridVolume *vol = new GridVolume(color_t(sa), color_t(ss), color_t(le), g,
									point3d_t(min[0], min[1], min[2]),
									point3d_t(max[0], max[1], max[2]),
									attSc, densityFile);

	return vol;
}

extern "C"
{
	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("GridVolume", GridVolume::factory);
	}
}

__END_YAFRAY
