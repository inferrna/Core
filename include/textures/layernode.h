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

#ifndef Y_LAYERNODE_H
#define Y_LAYERNODE_H

#include <core_api/shader.h>
#include <core_api/texture.h>
#include <core_api/environment.h>

__BEGIN_YAFRAY

// texture flag
#define TXF_RGBTOINT    1
#define TXF_STENCIL     2
#define TXF_NEGATIVE    4
#define TXF_ALPHAMIX    8

class layerNode_t: public shaderNode_t
{
	public:
		layerNode_t(unsigned tflag, CFLOAT col_fac, CFLOAT var_fac, CFLOAT def_val, colorA_t def_col, mix_modes mmod);
		virtual void eval(nodeStack_t &stack, const renderState_t &state, const surfacePoint_t &sp)const;
		virtual void eval(nodeStack_t &stack, const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wi)const;
		virtual void evalDerivative(nodeStack_t &stack, const renderState_t &state, const surfacePoint_t &sp)const;
		virtual bool isViewDependant() const;
		virtual bool configInputs(const paraMap_t &params, const nodeFinder_t &find);
		//virtual void getDerivative(const surfacePoint_t &sp, float &du, float &dv)const;
		virtual bool getDependencies(std::vector<const shaderNode_t*> &dep) const;
		static shaderNode_t* factory(const paraMap_t &params,renderEnvironment_t &render);
	protected:
		const shaderNode_t *input, *upperLayer;
		unsigned int texflag;
		CFLOAT colfac;
		CFLOAT valfac;
		CFLOAT default_val, upper_val;
		colorA_t default_col, upper_col;
		mix_modes mode;
		bool do_color, do_scalar, color_input, use_alpha;
};



__END_YAFRAY

#endif // Y_LAYERNODE_H
