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

#ifndef Y_SAMPLING_H
#define Y_SAMPLING_H

#include<yafray_config.h>

#include "imagesplitter.h"
#include <vector>

__BEGIN_YAFRAY


class sampler_t
{
	public:
	sampler_t(int width, int height);
	/*! add dimension to the sampler
		\return reached dimension, needed to request the samples */
	int addDimension(int max_samples);
	/*! set the image area to generate samples for; remember that each thread
		has its own sampler instance, which gets assigned image areas this way
		until the whole image has been rendered.
	*/
	bool setArea(const renderArea_t &a);
	//! get the next value in the sequence for dimension
	float getNext(unsigned int dimension);
//	bool getSample(sample_t &s);
	protected:
	std::vector<sequence_t> dim; //!< strore the sequences for each dimension
}



__END_YAFRAY

#endif // Y_SAMPLING_H
