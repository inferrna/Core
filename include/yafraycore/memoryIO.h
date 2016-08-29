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

#ifndef Y_MEMORYIO_H
#define Y_MEMORYIO_H

#include <core_api/color.h>
#include <utilities/buffer.h>
#include <core_api/output.h>

__BEGIN_YAFRAY

class YAFRAYCORE_EXPORT memoryIO_t : public colorOutput_t
{
	public:
		memoryIO_t(int resx, int resy, float* iMem);
		virtual bool putPixel(int x, int y, const float *c, bool alpha = true, bool depth = false, float z = 0.f);
		void flush();
		virtual void flushArea(int x0, int y0, int x1, int y1) {}; // no tiled file format used...yet
		virtual ~memoryIO_t();
	protected:
		int sizex, sizey;
		float* imageMem;
};


__END_YAFRAY

#endif // Y_MEMORYIO_H
