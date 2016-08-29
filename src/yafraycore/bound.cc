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

#include <core_api/bound.h>
#include <iostream>

__BEGIN_YAFRAY

bound_t::bound_t(const bound_t &r,const bound_t &l)
{
	PFLOAT minx=std::min(r.a.x,l.a.x);
	PFLOAT miny=std::min(r.a.y,l.a.y);
	PFLOAT minz=std::min(r.a.z,l.a.z);
	PFLOAT maxx=std::max(r.g.x,l.g.x);
	PFLOAT maxy=std::max(r.g.y,l.g.y);
	PFLOAT maxz=std::max(r.g.z,l.g.z);
	a.set(minx,miny,minz);
	g.set(maxx,maxy,maxz);
}

GFLOAT bound_t::vol() const
{
	GFLOAT ret=(g.y-a.y)*(g.x-a.x)*(g.z-a.z);

	return ret;
}

__END_YAFRAY
