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

#ifndef Y_QTAPI_H
#define Y_QTAPI_H

#ifdef BUILDING_QTPLUGIN
  #define YAFARAY_QT_EXPORT YF_EXPORT
#else
  #define YAFARAY_QT_EXPORT YF_IMPORT
#endif
#include <interface/yafrayinterface.h>
#include <string>

namespace yafaray
{
	class yafrayInterface_t;
}

struct YAFARAY_QT_EXPORT Settings {
	bool autoSave;
	bool autoSaveAlpha;
	bool closeAfterFinish;
	std::string fileName;
};

extern "C"
{
	YAFARAY_QT_EXPORT void initGui();
	YAFARAY_QT_EXPORT int createRenderWidget(yafaray::yafrayInterface_t *interf, int xsize, int ysize, int bStartX, int bStartY, Settings settings);
}

#endif // Y_QTAPI_H

