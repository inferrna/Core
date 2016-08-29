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

#ifndef __YAFSYSTEM_H
#define __YAFSYSTEM_H

#include <yafray_config.h>

#ifdef WIN32
	#include <io.h>
	#include <windows.h>
#endif

#include <list>
#include <string>

__BEGIN_YAFRAY

class YAFRAYCORE_EXPORT sharedlibrary_t 
{
	public:
	  sharedlibrary_t();
	  sharedlibrary_t(const std::string &library);
	  sharedlibrary_t(const sharedlibrary_t &src);
	  ~sharedlibrary_t();

	  bool isOpen();
	  void* getSymbol(const char *name);

	protected:

	  void open(const std::string &library);
	  void close();
  	void addReference() { (*refcount)++; };
  	void removeReference() { (*refcount)--; };
	bool isUsed()const {return ((*refcount)>0);}; 


	int *refcount;
#ifdef WIN32
	HINSTANCE handle;
#else
	void *handle;
#endif
};

YAFRAYCORE_EXPORT const std::list<std::string> & listDir(const std::string &dir);


__END_YAFRAY

#endif
