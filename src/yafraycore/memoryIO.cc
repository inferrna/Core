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

#include <core_api/color.h>
#include <utilities/buffer.h>
#include <core_api/output.h>
#include <yafraycore/memoryIO.h>
#include <cstdlib>

__BEGIN_YAFRAY


memoryIO_t::memoryIO_t ( int resx, int resy, float* iMem )
{
	sizex = resx;
	sizey = resy;
	imageMem = iMem; // iMem must be a valid pointer to memory of the size: sizex * sizey * 4 * sizeof(float)
}

// Depth channel support?
bool memoryIO_t::putPixel ( int x, int y, const float *c, bool alpha, bool depth, float z )
{
	for (int i = 0; i < 4; ++i)
	{
		if(!alpha && i == 3) imageMem[(x + sizex * y) * 4 + i] = 1.f;
		else imageMem[(x + sizex * y) * 4 + i] = c[i];
	}
	return true;
}

void memoryIO_t::flush() { }

memoryIO_t::~memoryIO_t() { }


__END_YAFRAY

