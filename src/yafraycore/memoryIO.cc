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
bool memoryIO_t::putPixel(int numView, int x, int y, const renderPasses_t *renderPasses, const std::vector<colorA_t> &colExtPasses, bool alpha)
{
	imageMem[(x + sizex * y) * 4 + 0] = colExtPasses.at(0).R;
	imageMem[(x + sizex * y) * 4 + 0] = colExtPasses.at(0).G;
	imageMem[(x + sizex * y) * 4 + 0] = colExtPasses.at(0).B;
	if(!alpha) imageMem[(x + sizex * y) * 4 + 3] = 1.f;

	return true;
}

void memoryIO_t::flush(int numView, const renderPasses_t *renderPasses) { }

memoryIO_t::~memoryIO_t() { }


__END_YAFRAY

