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

#ifndef Y_PHOTONMAP_H
#define Y_PHOTONMAP_H

#include <yafray_config.h>
// povman: photonkdtree.h is a more advanced code added for SSS translucent material..
// contain an expanded code from pkdtree.h
#include "photonkdtree.h"
// end
#include "pkdtree.h"
#include <core_api/color.h>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

__BEGIN_YAFRAY
#define c255Ratio 81.16902097686662123083
#define c256Ratio 40.74366543152520595687


#define cInv255Ratio 0.01231997119054820878
#define cInv256Ratio 0.02454369260617025968


class dirConverter_t
{
	public:
		dirConverter_t();

		vector3d_t convert(unsigned char theta,unsigned char phi)
		{
			return vector3d_t(sintheta[theta]*cosphi[phi],
								sintheta[theta]*sinphi[phi],
								costheta[theta]);
		}
		std::pair<unsigned char,unsigned char> convert(const vector3d_t &dir)
		{
			int t=(int)(fAcos(dir.z)*c255Ratio);
			int p=(int)(atan2(dir.y,dir.x)*c256Ratio);
			if(t>254) t=254;
			else if(t<0) t=0;
			if(p>255) p=255;
			else if(p<0) p+=256;
			return std::pair<unsigned char,unsigned char>(t,p);
		}

	protected:
		float cosphi[256];
		float sinphi[256];
		float costheta[255];
		float sintheta[255];
};

extern YAFRAYCORE_EXPORT dirConverter_t dirconverter;

class photon_t
{
	public:
		photon_t() {/*theta=255;*/};
		photon_t(const vector3d_t &d,const point3d_t &p, const color_t &col)
		{
#ifdef _SMALL_PHOTONS
			direction(d);
#else
			dir=d;
#endif
			pos=p;
			c=col;
		};
//		photon_t(const runningPhoton_t &p)
//		{
//			pos=p.pos;
//			c=p.c;
//			vector3d_t dir=p.lastpos-p.pos;
//			dir.normalize();
//			direction(dir);
//		};
		const point3d_t & position()const {return pos;};
		const color_t color()const {return c;};
		void color(const color_t &col) {c=col;};
		vector3d_t direction()const
		{
#ifdef _SMALL_PHOTONS
			if(theta==255) return vector3d_t(0,0,0);
			else return dirconverter.convert(theta,phi);
#else
			return (vector3d_t)dir;
#endif
		};
		void direction(const vector3d_t &d)
		{
#ifdef _SMALL_PHOTONS
			if(dir.null()) theta=255;
			else
			{
				std::pair<unsigned char,unsigned char> cd=dirconverter.convert(d);
				theta=cd.first;
				phi=cd.second;
			}
#else
			dir=d;
#endif
		}

		point3d_t pos;
#ifdef _SMALL_PHOTONS
		rgbe_t c;
		unsigned char theta,phi;

		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & BOOST_SERIALIZATION_NVP(pos);
			ar & BOOST_SERIALIZATION_NVP(c);
			ar & BOOST_SERIALIZATION_NVP(theta);
			ar & BOOST_SERIALIZATION_NVP(phi);
		}
#else
		color_t c;
		normal_t dir;

		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & BOOST_SERIALIZATION_NVP(pos);
			ar & BOOST_SERIALIZATION_NVP(c);
			ar & BOOST_SERIALIZATION_NVP(dir);
		}
#endif
        vector3d_t hitNormal;

		// for photon inner sss object
		point3d_t sourcePos;
		float	  sourceDepth;
};

struct radData_t
{
	radData_t(point3d_t &p, vector3d_t n): pos(p), normal(n), use(true) {}
	point3d_t pos;
	vector3d_t normal;
	color_t refl;
	color_t transm;
	mutable bool use;
};

struct foundPhoton_t
{
	foundPhoton_t(){};
	foundPhoton_t(const photon_t *p, float d): photon(p), distSquare(d){}
	bool operator<(const foundPhoton_t &p2) const { return distSquare < p2.distSquare; }
	const photon_t *photon;
	float distSquare;
	//temp!!
	float dis;
};

class YAFRAYCORE_EXPORT photonMap_t
{
public:
    // povman: for SSS, expand constructor and destructor for add phtree

    //photonMap_t(): paths(0), updated(false), searchRadius(1.), tree(0){ }
    photonMap_t(): paths(0), updated(false), searchRadius(1.), tree(0), phTree(0){ }

    // expand destructor are in photon.cc
    //~photonMap_t(){ if(tree) delete tree; }
    ~photonMap_t();

    // end
    void setNumPaths(int n){ paths=n; }
    int nPaths() const{ return paths; }
    int nPhotons() const{ return photons.size(); }
    void pushPhoton(photon_t &p) { photons.push_back(p); updated=false; }
    void swapVector(std::vector<photon_t> &vec) { photons.swap(vec); updated=false; }
    void updateTree();
    // povman: for SSS, expand function clear for add phTree
    //void clear(){ photons.clear(); delete tree; tree=0; updated=false; }
    void clear(){ photons.clear(); delete tree; tree=0; delete phTree; phTree=0; updated=false; }
    // end
    bool ready() const { return updated; }
    //void gather(const point3d_t &P, std::vector< foundPhoton_t > &found, unsigned int K, PFLOAT &sqRadius) const;
    int gather(const point3d_t &P, foundPhoton_t *found, unsigned int K, PFLOAT &sqRadius) const;
    const photon_t* findNearest(const point3d_t &P, const vector3d_t &n, PFLOAT dist) const;

    // add SSS code
    void getAllPhotons(const point3d_t& woP, std::vector<const photon_t*>& sssPhotons);
    int numberOfPhotonInDisc(const point3d_t &p, PFLOAT scale, PFLOAT dist) const;
    // test
    void updatePhTree();

    // end
protected:
    std::vector<photon_t> photons;
    int paths; //!< amount of photon paths that have been traced for generating the map
    bool updated;
    PFLOAT searchRadius;
    kdtree::pointKdTree<photon_t> *tree;
    // povman: add photonKdTree for use with SSS maps
    kdtree::photonKdTree<photon_t> *phTree;
};

// photon "processes" for lookup

struct photonGather_t
{
	photonGather_t(u_int32 mp, const point3d_t &p);
	void operator()(const photon_t *photon, float dist2, float &maxDistSquared) const;
	const point3d_t &p;
	foundPhoton_t *photons;
	u_int32 nLookup;
	mutable u_int32 foundPhotons;
};

struct nearestPhoton_t
{
	nearestPhoton_t(const point3d_t &pos, const vector3d_t &norm): p(pos), n(norm), nearest(nullptr) {}
	void operator()(const photon_t *photon, float dist2, float &maxDistSquared) const
	{
		if ( photon->direction() * n > 0.f) { nearest = photon; maxDistSquared = dist2; }
	}
	const point3d_t p; //wth do i need this for actually??
	const vector3d_t n;
	mutable const photon_t *nearest;
};

/*! "eliminates" photons within lookup radius (sets use=false) */
struct eliminatePhoton_t
{
	eliminatePhoton_t(const vector3d_t &norm): n(norm) {}
	void operator()(const radData_t *rpoint, float dist2, float &maxDistSquared) const
	{
		if ( rpoint->normal * n > 0.f) { rpoint->use = false; }
	}
	const vector3d_t n;
};


YAFRAYCORE_EXPORT bool photonMapLoad(photonMap_t * map, const std::string &filename, bool debugXMLformat = false);

YAFRAYCORE_EXPORT bool photonMapSave(const photonMap_t * map, const std::string &filename, bool debugXMLformat = false);


__END_YAFRAY

#endif // Y_PHOTONMAP_H
