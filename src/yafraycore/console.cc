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

#include <yafraycore/monitor.h>
#include <utilities/math_utils.h>
#include <iostream>
#include <string>
#include <iomanip>
#ifdef _MSC_VER
#include <algorithm>
#endif

__BEGIN_YAFRAY

#define printBar(progEmpty, progFull, per) \
std::cout << "\r"; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(Green); \
std::cout << "Progress: "; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(Red, true); \
std::cout << "["; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(Green, true); \
std::cout << std::string(progFull, '#') << std::string(progEmpty, ' '); \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(Red, true); \
std::cout << "] "; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(); \
std::cout << "("; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(Yellow, true); \
std::cout << per << "%"; \
if(yafLog.getConsoleLogColorsEnabled()) std::cout << setColor(); \
std::cout << ")" << std::flush

ConsoleProgressBar_t::ConsoleProgressBar_t(int cwidth): width(cwidth), nSteps(0), doneSteps(0)
{
	totalBarLen = width - 22;
}

void ConsoleProgressBar_t::init(int totalSteps)
{
	nSteps=totalSteps;
	doneSteps = 0;
	lastBarLen = 0;
	printBar(totalBarLen, 0, 0);
}

void ConsoleProgressBar_t::update(int steps)
{
	doneSteps += steps;
	float progress = (float) std::min(doneSteps, nSteps) / (float) nSteps;
	int barLen = std::min(totalBarLen, (int)(totalBarLen*progress));
	if(!(barLen >= 0)) barLen = 0;
	if(barLen > lastBarLen)
	{
		printBar(totalBarLen-barLen, barLen, (int) (100 * progress));
	}
	lastBarLen = barLen;
}

void ConsoleProgressBar_t::done()
{
	printBar(0, totalBarLen, 100) << yendl;
}

float ConsoleProgressBar_t::getPercent() const
{
	float progress = 100.f * RoundFloatPrecision((float) std::min(doneSteps, nSteps) / (float) nSteps, 0.01);
	return progress;
}

__END_YAFRAY
