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

#ifndef Y_BACKGROUND_H
#define Y_BACKGROUND_H

#include <yafray_config.h>

#include "color.h"
#include "ray.h"

__BEGIN_YAFRAY

struct renderState_t;
class light_t;

class YAFRAYCORE_EXPORT background_t
{
	public:
		//! get the background color for a given ray
		virtual color_t operator() (const ray_t &ray, renderState_t &state, bool filtered=false) const=0;
		virtual color_t eval(const ray_t &ray, bool filtered=false) const=0;
		/*! get the light source representing background lighting.
			\return the light source that reproduces background lighting, or NULL if background
					shall only be sampled from BSDFs
		*/
		virtual ~background_t() {};
};

__END_YAFRAY

#endif // Y_BACKGROUND_H
