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

#ifndef Y_MATHUTIL_H
#define Y_MATHUTIL_H

#include <cmath>

/* povman..
move repeat code from kdtree.cc and ray_kdtree.cc
and add #include <math_utils.h> to kdtree headers files
*/
#if (defined(_M_IX86) || defined(i386) || defined(_X86_))
#define Y_FAST_INT 1
#else
#define Y_FAST_INT 0
#endif

#define _doublemagicroundeps    (0.5 - 1.4e-11)
//almost .5f = .5f - 1e^(number of exp bit)

#define _doublemagic            double (6755399441055744.0)
//2^52 * 1.5,  uses limited precision to floor

inline int Y_Round2Int(double val) {
#if Y_FAST_INT > 0
    union { double f; int i[2]; } u;
    u.f = val + _doublemagic;
    return u.i[0];
#else
    return int(val);
#endif
}

inline int Y_Float2Int(double val) {
#if Y_FAST_INT > 0
    return (val<0) ? Y_Round2Int(val + _doublemagicroundeps) :
        Y_Round2Int(val - _doublemagicroundeps);
#else
    return (int)val;
#endif
}
//----------------------------------------------------

 /* add Y_FAST_INT flag to old code. FAST_INT is already undefined.. */

inline int Round2Int(double val) {
#if Y_FAST_INT > 0
    val = val + _doublemagic;
        return ((long*)&val)[0];
#else
//  #warning "using slow rounding"
        return int (val+_doublemagicroundeps);
#endif
}

inline int Float2Int(double val) {
#if Y_FAST_INT > 0
    return (val<0) ?  Round2Int(val+_doublemagicroundeps) : Round2Int(val-_doublemagicroundeps);
#else
//  #warning "using slow rounding"
    return (int)val;
#endif
}

inline int Floor2Int(double val) {
#if Y_FAST_INT > 0
    return Round2Int(val - _doublemagicroundeps);
#else
//  #warning "using slow rounding"
    return (int)std::floor(val);
#endif
}

inline int Ceil2Int(double val) {
#if Y_FAST_INT > 0
    return Round2Int(val + _doublemagicroundeps);
#else
//  #warning "using slow rounding"
    return (int)std::ceil(val);
#endif
}
//povman: add in range inline fuction
inline float inRange(float up, float down, float val){
    if (val > up) val = up;
    if (val < down) val = down;
    return val;
}

#endif // Y_MATHUTIL_H
