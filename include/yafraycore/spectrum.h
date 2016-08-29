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

#ifndef Y_SPECTRUM_H
#define Y_SPECTRUM_H

#include <core_api/color.h>

__BEGIN_YAFRAY

YAFRAYCORE_EXPORT void wl2rgb_fromCIE(float wl, color_t &col);
//YAFRAYCORE_EXPORT void approxSpectrumRGB(float wl, color_t &col);
//YAFRAYCORE_EXPORT void fakeSpectrum(float p, color_t &col);
YAFRAYCORE_EXPORT void CauchyCoefficients(float IOR, float disp_pw, float &CauchyA, float &CauchyB);
YAFRAYCORE_EXPORT float getIORcolor(float w, float CauchyA, float CauchyB, color_t &col);
YAFRAYCORE_EXPORT color_t wl2XYZ(float wl);

static inline float getIOR(float w, float CauchyA, float CauchyB)
{
	float wl = 300.0*w + 400.0;
	return CauchyA + CauchyB/(wl*wl);
}

static inline void wl2rgb(float w, color_t &wl_col)
{
	float wl = 300.0*w + 400.0;
	wl2rgb_fromCIE(wl, wl_col);
	wl_col *= 2.214032659670777114f;
}

__END_YAFRAY

#endif // Y_SPECTRUM_H
