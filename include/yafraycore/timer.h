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

#ifndef Y_TIMER_H
#define Y_TIMER_H

#include <yafray_config.h>

#include <time.h>
#include <string>
#include <map>

#ifndef _WIN32
extern "C" struct timeval;
#include <sys/time.h>
#endif

__BEGIN_YAFRAY

class YAFRAYCORE_EXPORT timer_t
{
	public:
		bool addEvent(const std::string &name);
		bool start(const std::string &name);
		bool stop(const std::string &name);
		bool reset(const std::string &name);
		double getTime(const std::string &name);
		
		static void splitTime(double t, double *secs, int *mins=0, int *hours=0, int *days=0);
	
	protected:
		bool includes(const std::string &label)const;
		
		struct tdata_t
		{
			tdata_t():started(false), stopped(false) {};
			clock_t start, finish;
			#ifndef WIN32
			timeval tvs, tvf;
			#endif
			bool started, stopped;
		};
		std::map<std::string, tdata_t> events;
};

// global timer object, defined in timer.cc
extern YAFRAYCORE_EXPORT timer_t gTimer;

__END_YAFRAY

#endif
