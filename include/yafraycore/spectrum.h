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

YAFRAYCORE_EXPORT void wl2rgb_fromCIE(CFLOAT wl, color_t &col);
//YAFRAYCORE_EXPORT void approxSpectrumRGB(CFLOAT wl, color_t &col);
//YAFRAYCORE_EXPORT void fakeSpectrum(CFLOAT p, color_t &col);
YAFRAYCORE_EXPORT void CauchyCoefficients(PFLOAT IOR, PFLOAT disp_pw, PFLOAT &CauchyA, PFLOAT &CauchyB);
YAFRAYCORE_EXPORT PFLOAT getIORcolor(PFLOAT w, PFLOAT CauchyA, PFLOAT CauchyB, color_t &col);
YAFRAYCORE_EXPORT color_t wl2XYZ(CFLOAT wl);

static inline PFLOAT getIOR(PFLOAT w, PFLOAT CauchyA, PFLOAT CauchyB)
{
	PFLOAT wl = 300.0*w + 400.0;
	return CauchyA + CauchyB/(wl*wl);
}

static inline void wl2rgb(float w, color_t &wl_col)
{
	PFLOAT wl = 300.0*w + 400.0;
	wl2rgb_fromCIE(wl, wl_col);
	wl_col *= 2.214032659670777114f;
}

__END_YAFRAY

#endif // Y_SPECTRUM_H
