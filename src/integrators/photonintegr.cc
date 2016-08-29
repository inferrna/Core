/****************************************************************************
 *      photonintegr.cc: integrator for photon mapping and final gather
 *      This is part of the yafaray package
 *      Copyright (C) 2006  Mathias Wein (Lynx)
 *		Copyright (C) 2009  Rodrigo Placencia (DarkTide)
 *
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

#include <integrators/photonintegr.h>
#include <utilities/mcqmc.h>
#include <yafraycore/scr_halton.h>

#include <sstream>
#include <iomanip>

__BEGIN_YAFRAY

struct preGatherData_t
{
	preGatherData_t(photonMap_t *dm): diffuseMap(dm), fetched(0) {}
	photonMap_t *diffuseMap;

	std::vector<radData_t> rad_points;
	std::vector<photon_t> radianceVec;
	progressBar_t *pbar;
	volatile int fetched;
	yafthreads::mutex_t mutex;
};

class preGatherWorker_t: public yafthreads::thread_t
{
	public:
		preGatherWorker_t(preGatherData_t *dat, float dsRad, int search):
			gdata(dat), dsRadius_2(dsRad*dsRad), nSearch(search) {};
		virtual void body();
	protected:
		preGatherData_t *gdata;
		float dsRadius_2;
		int nSearch;
};

void preGatherWorker_t::body()
{
	unsigned int start, endWork, total;

	gdata->mutex.lock();
	start = gdata->fetched;
	total = gdata->rad_points.size();
	endWork = gdata->fetched = std::min(total, start + 32);
	gdata->mutex.unlock();

	foundPhoton_t *gathered = new foundPhoton_t[nSearch];

	float radius = 0.f;
	float iScale = 1.f / ((float)gdata->diffuseMap->nPaths() * M_PI);
	float scale = 0.f;

	while(start < total)
	{
		for(unsigned int n=start; n<endWork; ++n)
		{
			radius = dsRadius_2;//actually the square radius...
			int nGathered = gdata->diffuseMap->gather(gdata->rad_points[n].pos, gathered, nSearch, radius);

			vector3d_t rnorm = gdata->rad_points[n].normal;

			color_t sum(0.0);

			if(nGathered > 0)
			{
				scale = iScale / radius;

				for(int i=0; i<nGathered; ++i)
				{
					vector3d_t pdir = gathered[i].photon->direction();

					if( rnorm * pdir > 0.f ) sum += gdata->rad_points[n].refl * scale * gathered[i].photon->color();
					else sum += gdata->rad_points[n].transm * scale * gathered[i].photon->color();
				}
			}

			gdata->radianceVec[n] = photon_t(rnorm, gdata->rad_points[n].pos, sum);
		}
		gdata->mutx.lock();
		start = gdata->fetched;
		endWork = gdata->fetched = std::min(total, start + 32);
		gdata->pbar->update(32);
		gdata->mutx.unlock();
	}
	delete[] gathered;
}

photonIntegrator_t::photonIntegrator_t(unsigned int dPhotons, unsigned int cPhotons, bool transpShad, int shadowDepth, float dsRad, float cRad)
{
	usePhotonCaustics = true;
	usePhotonDiffuse = true;
	type = SURFACE;
	trShad = transpShad;
	finalGather = true;
	nDiffusePhotons = dPhotons;
	nCausPhotons = cPhotons;
	sDepth = shadowDepth;
	dsRadius = dsRad;
	causRadius = cRad;
	rDepth = 6;
	maxBounces = 5;
	integratorName = "PhotonMap";
	integratorShortName = "PM";
}

photonIntegrator_t::~photonIntegrator_t()
{
	destorySSSMaps();
}


void photonIntegrator_t::causticWorker(photonMap_t * causticMap, int threadID, const scene_t *scene, unsigned int nCausPhotons, const pdf1D_t *lightPowerD, int numCLights, const std::string &integratorName, const std::vector<light_t *> &tmplights, int causDepth, progressBar_t *pb, int pbStep, unsigned int &totalPhotonsShot, int maxBounces)
{
	ray_t ray;
	float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
	color_t pcol;

	//shoot photons
	bool done=false;
	unsigned int curr=0;
	
	surfacePoint_t sp;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	state.cam = scene->getCamera();

	float fNumLights = (float)numCLights;
	unsigned int nCausPhotons_thread = 1 + ( (nCausPhotons - 1) / scene->getNumThreadsPhotons() );

	std::vector<photon_t> localCausticPhotons;
	localCausticPhotons.clear();
	localCausticPhotons.reserve(nCausPhotons_thread);

	if(!set.str().empty()) set << "+";

	set << "DiffPhotons [" << nDiffusePhotons << "]+CausPhotons[" << nCausPhotons << "]";

	if(finalGather)
	{
		set << "+FG[" << nPaths << ", " << gatherBounces << "]";
	}

	settings = set.str();

	ray_t ray;
	float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
	int numCLights = 0;
	int numDLights = 0;
	float fNumLights = 0.f;
	float *energies = NULL;
	color_t pcol;

		state.chromatic = true;
		state.wavelength = scrHalton(5,haltoncurr);

		s1 = RI_vdC(haltoncurr);
		s2 = scrHalton(2, haltoncurr);
		s3 = scrHalton(3, haltoncurr);
		s4 = scrHalton(4, haltoncurr);

		sL = float(haltoncurr) * invCaustPhotons;
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		
		if(lightNum >= numCLights)
		{
			causticMap->mutx.lock();
			Y_ERROR << integratorName << ": lightPDF sample error! " << sL << "/" << lightNum << yendl;
			causticMap->mutx.unlock();
			return;
		}
	}

	fNumLights = (float)numDLights;
	energies = new float[numDLights];

		pcol = tmplights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
		ray.tmin = scene->rayMinDist;
		ray.tmax = -1.0;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
		if(pcol.isBlack())
		{
			++curr;
			done = (curr >= nCausPhotons_thread);
			continue;
		}
		int nBounces=0;
		bool causticPhoton = false;
		bool directPhoton = true;
		const material_t *material = nullptr;
		BSDF_t bsdfs;

		while( scene->intersect(ray, sp) )
		{
			if(std::isnan(pcol.R) || std::isnan(pcol.G) || std::isnan(pcol.B))
			{
				causticMap->mutx.lock();
				Y_WARNING << integratorName << ": NaN  on photon color for light" << lightNum + 1 << "." << yendl;
				causticMap->mutx.unlock();
				continue;
			}
			
			color_t transm(1.f);
			color_t vcol(0.f);
			const volumeHandler_t* vol = nullptr;
			
			if(material)
			{
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * -ray.dir < 0)))
				{
					if(vol->transmittance(state, ray, vcol)) transm = vcol;
				}
			}
			
			vector3d_t wi = -ray.dir, wo;
			material = sp.material;
			material->initBSDF(state, sp, bsdfs);

	lightPowerD = new pdf1D_t(energies, numDLights);

	Y_INFO << integratorName << ": Light(s) photon color testing for diffuse map:" << yendl;

	for(int i=0;i<numDLights;++i)
	{
		pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
		lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
		Y_INFO << integratorName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
	}

	delete[] energies;

	//shoot photons
	bool done=false;
	unsigned int curr=0;

	// for radiance map:
	preGatherData_t pgdat(&diffuseMap);

	surfacePoint_t sp;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	state.cam = scene->getCamera();
	progressBar_t *pb;
	int pbStep;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);

	Y_INFO << integratorName << ": Building diffuse photon map..." << yendl;

	pb->init(128);
	pbStep = std::max(1U, nDiffusePhotons / 128);
	pb->setTag("Building diffuse photon map...");
	//Pregather diffuse photons

	std::vector<photon_t> localDiffusePhotons;
	std::vector<radData_t> localRadPoints;

	localDiffusePhotons.clear();
	localDiffusePhotons.reserve(nDiffusePhotons_thread);
	localRadPoints.clear();
	
	float invDiffPhotons = 1.f / (float)nDiffusePhotons;

	while(!done)
	{
		unsigned int haltoncurr = curr + nDiffusePhotons_thread * threadID;

		s1 = RI_vdC(haltoncurr);
		s2 = scrHalton(2, haltoncurr);
		s3 = scrHalton(3, haltoncurr);
		s4 = scrHalton(4, haltoncurr);

		sL = float(haltoncurr) * invDiffPhotons;
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		if(lightNum >= numDLights)
		{
			diffuseMap->mutx.lock();
			Y_ERROR << integratorName << ": lightPDF sample error! " << sL << "/" << lightNum << yendl;
			diffuseMap->mutx.unlock();
			return;
		}

		pcol = tmplights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
		ray.tmin = scene->rayMinDist;
		ray.tmax = -1.0;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...

		if(pcol.isBlack())
		{
			++curr;
			done = (curr >= nDiffusePhotons_thread);
			continue;
		}

		int nBounces=0;
		bool causticPhoton = false;
		bool directPhoton = true;
		const material_t *material = nullptr;
		BSDF_t bsdfs;

		while( scene->intersect(ray, sp) )
		{
			if(std::isnan(pcol.R) || std::isnan(pcol.G) || std::isnan(pcol.B))
			{
				diffuseMap->mutx.lock();
				Y_WARNING << integratorName << ": NaN  on photon color for light" << lightNum + 1 << "." << yendl;
				diffuseMap->mutx.unlock();
				continue;
			}

			color_t transm(1.f);
			color_t vcol(0.f);
			const volumeHandler_t* vol = NULL;

			if(material)
			{
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * -ray.dir < 0)))
				{
					if(vol->transmittance(state, ray, vcol)) transm = vcol;
				}
			}

			vector3d_t wi = -ray.dir, wo;
			material = sp.material;
			material->initBSDF(state, sp, bsdfs);

			if(bsdfs & (BSDF_DIFFUSE))
			{
				//deposit photon on surface
				if(!causticPhoton)
				{
					photon_t np(wi, sp.P, pcol);
					localDiffusePhotons.push_back(np);
				}
				// create entry for radiance photon:
				// don't forget to choose subset only, face normal forward; geometric vs. smooth normal?
				if(finalGather && ourRandom() < 0.125 && !causticPhoton )
				{
					vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wi);
					radData_t rd(sp.P, N);
					rd.refl = material->getReflectivity(state, sp, BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECT);
					rd.transm = material->getReflectivity(state, sp, BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_TRANSMIT);
					localRadPoints.push_back(rd);
				}
			}
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces) break;
			// scatter photon
			int d5 = 3*nBounces + 5;

			s5 = scrHalton(d5, curr);
			s6 = scrHalton(d5+1, curr);
			s7 = scrHalton(d5+2, curr);

			pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

			bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
			if(!scattered) break; //photon was absorbed.

			pcol = sample.color;

			causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
							((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
			directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;

			ray.from = sp.P;
			ray.dir = wo;
			ray.tmin = scene->rayMinDist;
			ray.tmax = -1.0;
			++nBounces;
		}
		++curr;
		if(curr % pbStep == 0)
		{
			pb->mutx.lock();
			pb->update();
			pb->mutx.unlock();
			if(scene->getSignals() & Y_SIG_ABORT) { return; }
		}
		done = (curr >= nDiffusePhotons_thread);
	}
	diffuseMap->mutx.lock();
	diffuseMap->appendVector(localDiffusePhotons, curr);
	totalPhotonsShot += curr;
	diffuseMap->mutx.unlock();
	
	pgdat.mutx.lock();
	pgdat.rad_points.insert(std::end(pgdat.rad_points), std::begin(localRadPoints), std::end(localRadPoints));
	pgdat.mutx.unlock();
}

void photonIntegrator_t::photonMapKdTreeWorker(photonMap_t * photonMap)
{
	photonMap->updateTree();
}

bool photonIntegrator_t::preprocess()
{
	progressBar_t *pb;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);

	lookupRad = 4*dsRadius*dsRadius;
		
	std::stringstream set;
	gTimer.addEvent("prepass");
	gTimer.start("prepass");

	Y_INFO << integratorName << ": Starting preprocess..." << yendl;

	set << "Photon Mapping  ";

	if(trShad)
	{
		set << "ShadowDepth=" << sDepth << "  ";
	}
	set << "RayDepth=" << rDepth << "  ";

	background = scene->getBackground();
	lights = scene->lights;
	std::vector<light_t*> tmplights;

	if(usePhotonCaustics)
	{
		set << "\nCaustic photons=" << nCausPhotons << " search=" << nCausSearch <<" radius=" << causRadius << " depth=" << causDepth << "  ";
	}

	if(usePhotonDiffuse)
	{
		set << "\nDiffuse photons=" << nDiffusePhotons << " search=" << nDiffuseSearch <<" radius=" << dsRadius << "  ";
	}

	if(finalGather)
	{
		set << " FG paths=" << nPaths << " bounces=" << gatherBounces << "  ";
	}
		
	if(photonMapProcessing == PHOTONS_LOAD)
	{
		bool causticMapFailedLoad = false;
		bool diffuseMapFailedLoad = false;
		bool fgRadianceMapFailedLoad = false;
		
		if(usePhotonCaustics)
		{
			pb->setTag("Loading caustic photon map from file...");
			std::string filename = session.getPathImageOutput() + "_caustic.photonmap";
			Y_INFO << integratorName << ": Loading caustic photon map from: " << filename << ". If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!" << yendl;
			if(photonMapLoad(session.causticMap, filename)) Y_VERBOSE << integratorName << ": Caustic map loaded." << yendl;
			else causticMapFailedLoad = true;
		}

		if(usePhotonDiffuse)
		{
			pb->setTag("Loading diffuse photon map from file...");
			std::string filename = session.getPathImageOutput() + "_diffuse.photonmap";
			Y_INFO << integratorName << ": Loading diffuse photon map from: " << filename << ". If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!"  << yendl;
			if(photonMapLoad(session.diffuseMap, filename)) Y_VERBOSE << integratorName << ": Diffuse map loaded." << yendl;
			else diffuseMapFailedLoad = true;
		}

		if(usePhotonDiffuse && finalGather)
		{
			pb->setTag("Loading FG radiance photon map from file...");
			std::string filename = session.getPathImageOutput() + "_fg_radiance.photonmap";
			Y_INFO << integratorName << ": Loading FG radiance photon map from: " << filename << ". If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!"  << yendl;
			if(photonMapLoad(session.radianceMap, filename)) Y_VERBOSE << integratorName << ": FG radiance map loaded." << yendl;
			else fgRadianceMapFailedLoad = true;
		}
		
		if(causticMapFailedLoad || diffuseMapFailedLoad || fgRadianceMapFailedLoad)
		{
			photonMapProcessing = PHOTONS_GENERATE_AND_SAVE;
			Y_WARNING << integratorName << ": photon maps loading failed, changing to Generate and Save mode." << yendl;
		}
	}

	if(photonMapProcessing == PHOTONS_REUSE)
	{
		if(usePhotonCaustics)
		{
			Y_INFO << integratorName << ": Reusing caustics photon map from memory. If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!" << yendl;
			if(session.causticMap->nPhotons() == 0)
			{
				Y_WARNING << integratorName << ": Caustic photon map enabled but empty, cannot be reused: changing to Generate mode." << yendl;
				photonMapProcessing = PHOTONS_GENERATE_ONLY;
			}
		}

		if(usePhotonDiffuse)
		{
			Y_INFO << integratorName << ": Reusing diffuse photon map from memory. If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!" << yendl;
			if(session.diffuseMap->nPhotons() == 0)
			{
				Y_WARNING << integratorName << ": Diffuse photon map enabled but empty, cannot be reused: changing to Generate mode." << yendl;
				photonMapProcessing = PHOTONS_GENERATE_ONLY;
			}
		}

		if(finalGather)
		{
			Y_INFO << integratorName << ": Reusing FG radiance photon map from memory. If it does not match the scene you could have crashes and/or incorrect renders, USE WITH CARE!" << yendl;
			if(session.radianceMap->nPhotons() == 0)
			{
				Y_WARNING << integratorName << ": FG radiance photon map enabled but empty, cannot be reused: changing to Generate mode." << yendl;
				photonMapProcessing = PHOTONS_GENERATE_ONLY;
			}
		}
	}

	if(photonMapProcessing == PHOTONS_LOAD)
	{
		set << " (loading photon maps from file)";
	}
	else if(photonMapProcessing == PHOTONS_REUSE)
	{
		set << " (reusing photon maps from memory)";
	}
	else if(photonMapProcessing == PHOTONS_GENERATE_AND_SAVE) set << " (saving photon maps to file)";

	if(photonMapProcessing == PHOTONS_LOAD || photonMapProcessing == PHOTONS_REUSE)
	{
		gTimer.stop("prepass");
		Y_INFO << integratorName << ": Photonmap building time: " << std::fixed << std::setprecision(1) << gTimer.getTime("prepass") << "s" << yendl;

		set << " [" << std::fixed << std::setprecision(1) << gTimer.getTime("prepass") << "s" << "]";

		yafLog.appendRenderSettings(set.str());
		
		for (std::string line; std::getline(set, line, '\n');) Y_VERBOSE << line << yendl;
		
		return true;
	}

	session.diffuseMap->clear();
	session.diffuseMap->setNumPaths(0);
	session.diffuseMap->reserveMemory(nDiffusePhotons);
	session.diffuseMap->setNumThreadsPKDtree(scene->getNumThreadsPhotons());

	session.causticMap->clear();
	session.causticMap->setNumPaths(0);
	session.causticMap->reserveMemory(nCausPhotons);
	session.causticMap->setNumThreadsPKDtree(scene->getNumThreadsPhotons());

	session.radianceMap->clear();
	session.radianceMap->setNumPaths(0);
	session.radianceMap->setNumThreadsPKDtree(scene->getNumThreadsPhotons());

	ray_t ray;
	float lightNumPdf, lightPdf;
	int numCLights = 0;
	int numDLights = 0;
	float fNumLights = 0.f;
	float *energies = nullptr;
	color_t pcol;

	//shoot photons
	unsigned int curr=0;
	// for radiance map:
	preGatherData_t pgdat(session.diffuseMap);
	
	surfacePoint_t sp;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	state.cam = scene->getCamera();
	int pbStep;

	tmplights.clear();
    // TODO: povman, investigate caustic lights..
	for(int i=0;i<(int)lights.size();++i)
	{
		if(lights[i]->shootsDiffuseP())
		{
			numDLights++;
			tmplights.push_back(lights[i]);
		}
	}
	
	if(numDLights == 0)
	{
		Y_WARNING << integratorName << ": No lights found that can shoot diffuse photons, disabling Diffuse photon processing" << yendl;
		enableDiffuse(false);
	}
	
	if( usePhotonDiffuse )
	{
		fNumLights = (float)numDLights;
		energies = new float[numDLights];

		for(int i=0;i<numDLights;++i) energies[i] = tmplights[i]->totalEnergy().energy();

		lightPowerD = new pdf1D_t(energies, numDLights);
		
		Y_VERBOSE << integratorName << ": Light(s) photon color testing for diffuse map:" << yendl;
		for(int i=0;i<numDLights;++i)
		{
			pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_VERBOSE << integratorName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
		}
		
		delete[] energies;
		
		//shoot photons
		curr=0;
		
		Y_INFO << integratorName << ": Building diffuse photon map..." << yendl;
		
		pb->init(128);
		pbStep = std::max(1U, nDiffusePhotons / 128);
		pb->setTag("Building diffuse photon map...");
		//Pregather diffuse photons

		int nThreads = scene->getNumThreadsPhotons();

		nDiffusePhotons = std::max((unsigned int) nThreads, (nDiffusePhotons / nThreads) * nThreads); //rounding the number of diffuse photons so it's a number divisible by the number of threads (distribute uniformly among the threads). At least 1 photon per thread
		
		Y_PARAMS << integratorName << ": Shooting "<<nDiffusePhotons<<" photons across " << nThreads << " threads (" << (nDiffusePhotons / nThreads) << " photons/thread)"<< yendl;
		
		if(nThreads >= 2)
		{
			std::vector<std::thread> threads;
			for(int i=0; i<nThreads; ++i) threads.push_back(std::thread(&photonIntegrator_t::diffuseWorker, this, session.diffuseMap, i, scene, nDiffusePhotons, lightPowerD, numDLights, std::ref(integratorName), tmplights, pb, pbStep, std::ref(curr), maxBounces, finalGather, std::ref(pgdat)));
			for(auto& t : threads) t.join();
		}
		else
		{
			bool done=false;

			float invDiffPhotons = 1.f / (float)nDiffusePhotons;
			float s1, s2, s3, s4, s5, s6, s7, sL;
			
			while(!done)
			{
				if(scene->getSignals() & Y_SIG_ABORT) {  pb->done(); if(!intpb) delete pb; return false; }

				s1 = RI_vdC(curr);
				s2 = scrHalton(2, curr);
				s3 = scrHalton(3, curr);
				s4 = scrHalton(4, curr);

				sL = float(curr) * invDiffPhotons;
				int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
				if(lightNum >= numDLights)
				{
					Y_ERROR << integratorName << ": lightPDF sample error! " << sL << "/" << lightNum << "... stopping now." << yendl;
					delete lightPowerD;
					return false;
				}

				pcol = tmplights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
				ray.tmin = scene->rayMinDist;
				ray.tmax = -1.0;
				pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
				
				if(pcol.isBlack())
				{
					++curr;
					done = (curr >= nDiffusePhotons);
					continue;
				}

				int nBounces=0;
				bool causticPhoton = false;
				bool directPhoton = true;
				const material_t *material = nullptr;
				BSDF_t bsdfs;

				while( scene->intersect(ray, sp) )
				{
					if(std::isnan(pcol.R) || std::isnan(pcol.G) || std::isnan(pcol.B))
					{
						Y_WARNING << integratorName << ": NaN  on photon color for light" << lightNum + 1 << "." << yendl;
						continue;
					}
					
					color_t transm(1.f);
					color_t vcol(0.f);
					const volumeHandler_t* vol = nullptr;
					
					if(material)
					{
						if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * -ray.dir < 0)))
						{
							if(vol->transmittance(state, ray, vcol)) transm = vcol;
						}
					}
					
					vector3d_t wi = -ray.dir, wo;
					material = sp.material;
					material->initBSDF(state, sp, bsdfs);
					
					if(bsdfs & (BSDF_DIFFUSE))
					{
						//deposit photon on surface
						if(!causticPhoton)
						{
							photon_t np(wi, sp.P, pcol);
							session.diffuseMap->pushPhoton(np);
							session.diffuseMap->setNumPaths(curr);
						}
						// create entry for radiance photon:
						// don't forget to choose subset only, face normal forward; geometric vs. smooth normal?
						if(finalGather && ourRandom() < 0.125 && !causticPhoton )
						{
							vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wi);
							radData_t rd(sp.P, N);
							rd.refl = material->getReflectivity(state, sp, BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECT);
							rd.transm = material->getReflectivity(state, sp, BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_TRANSMIT);
							pgdat.rad_points.push_back(rd);
						}
					}
					// need to break in the middle otherwise we scatter the photon and then discard it => redundant
					if(nBounces == maxBounces) break;
					// scatter photon
					int d5 = 3*nBounces + 5;

					s5 = scrHalton(d5, curr);
					s6 = scrHalton(d5+1, curr);
					s7 = scrHalton(d5+2, curr);
					
					pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

					bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
					if(!scattered) break; //photon was absorped.

					pcol = sample.color;

					causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
									((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
					directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;

					ray.from = sp.P;
					ray.dir = wo;
					ray.tmin = scene->rayMinDist;
					ray.tmax = -1.0;
					++nBounces;
				}
				++curr;
				if(curr % pbStep == 0) pb->update();
				done = (curr >= nDiffusePhotons);
			}
		}

		pb->done();
		pb->setTag("Diffuse photon map built.");
		Y_VERBOSE << integratorName << ": Diffuse photon map built." << yendl;
		Y_INFO << integratorName << ": Shot "<<curr<<" photons from " << numDLights << " light(s)" << yendl;

		delete lightPowerD;

		tmplights.clear();
		
		if(session.diffuseMap->nPhotons() < 50)
		{
			Y_ERROR << integratorName << ": Too few diffuse photons, stopping now." << yendl;
			return false;
		}

		Y_VERBOSE << integratorName << ": Stored diffuse photons: " << session.diffuseMap->nPhotons() << yendl;
	}
	else
	{
		Y_INFO << integratorName << ": Diffuse photon mapping disabled, skipping..." << yendl;
	}
	
	std::thread * diffuseMapBuildKdTree_thread = nullptr;
	
	if( usePhotonDiffuse && session.diffuseMap->nPhotons() > 0 && scene->getNumThreadsPhotons() >= 2)
	{
		Y_INFO << integratorName << ": Building diffuse photons kd-tree:" << yendl;
		pb->setTag("Building diffuse photons kd-tree...");

		diffuseMapBuildKdTree_thread = new std::thread(&photonIntegrator_t::photonMapKdTreeWorker, this, session.diffuseMap);
	}
	else

	if( usePhotonDiffuse && session.diffuseMap->nPhotons() > 0)
	{
		Y_INFO << integratorName << ": Building diffuse photons kd-tree:" << yendl;
		pb->setTag("Building diffuse photons kd-tree...");
		session.diffuseMap->updateTree();
		Y_VERBOSE << integratorName << ": Done." << yendl;
	}

	for(int i=0;i<(int)lights.size();++i)
	{
		if(lights[i]->shootsCausticP())
		{
			numCLights++;
			tmplights.push_back(lights[i]);
		}
	}
    // povman: if have.. proceed
	if(numCLights > 0)
	{
		done = false;
		curr=0;

		fNumLights = (float)numCLights;
		energies = new float[numCLights];
        // povman: sum of total energy..
		for(int i=0;i<numCLights;++i) energies[i] = tmplights[i]->totalEnergy().energy();

		lightPowerD = new pdf1D_t(energies, numCLights);

		Y_INFO << integratorName << ": Light(s) photon color testing for caustics map:" << yendl;
		for(int i=0;i<numCLights;++i)
		{
			pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_VERBOSE << integratorName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
		}

		delete[] energies;

		Y_INFO << integratorName << ": Building caustics photon map..." << yendl;
		pb->init(128);
		pbStep = std::max(1U, nCausPhotons / 128);
		pb->setTag("Building caustics photon map...");
		//Pregather caustic photons

		float invCaustPhotons = 1.f / (float)nCausPhotons;

		while(!done)
		{
			if(scene->getSignals() & Y_SIG_ABORT) { pb->done(); if(!intpb) delete pb; return false; }
			state.chromatic = true;
			state.wavelength = scrHalton(5,curr);


			sL = float(curr) * invCaustPhotons;
			int lightNum = lightPowerD->DSample(sL, &lightNumPdf);

			if(lightNum >= numCLights)
			{
				Y_ERROR << integratorName <<": lightPDF sample error! "<< sL <<"/"<< lightNum <<"... stopping now."<< yendl;
				delete lightPowerD;
				return false;
			}

			while(!done)
			{
				if(scene->getSignals() & Y_SIG_ABORT) { pb->done(); if(!intpb) delete pb; return false; }
				state.chromatic = true;
				state.wavelength = scrHalton(5,curr);

				s1 = RI_vdC(curr);
				s2 = scrHalton(2, curr);
				s3 = scrHalton(3, curr);
				s4 = scrHalton(4, curr);

				sL = float(curr) * invCaustPhotons;
				int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
				
				if(lightNum >= numCLights)
				{
					Y_WARNING << integratorName << ": NaN on photon color for light" << lightNum + 1 << "." << yendl;
					continue;
				}

				color_t transm(1.f);
				color_t vcol(0.f);
				const volumeHandler_t* vol = NULL;

				if(material)
				{
					if(std::isnan(pcol.R) || std::isnan(pcol.G) || std::isnan(pcol.B))
					{
						Y_WARNING << integratorName << ": NaN  on photon color for light" << lightNum + 1 << "." << yendl;
						continue;
					}
				}

				vector3d_t wi = -ray.dir, wo;
				material = sp.material;
				material->initBSDF(state, sp, bsdfs);

					if(bsdfs & BSDF_DIFFUSE)
					{
						if(causticPhoton)
						{
							photon_t np(wi, sp.P, pcol);
							session.causticMap->pushPhoton(np);
							session.causticMap->setNumPaths(curr);
						}
					}
				}

				// need to break in the middle otherwise we scatter the photon and then discard it => redundant
				if(nBounces == maxBounces) break;
				// scatter photon
				int d5 = 3*nBounces + 5;

					s5 = scrHalton(d5, curr);
					s6 = scrHalton(d5+1, curr);
					s7 = scrHalton(d5+2, curr);

					pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

				bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
				if(!scattered) break; //photon was absorbed.

					pcol = sample.color;

				causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
								((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
				directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;

				if(state.chromatic && (sample.sampledFlags & BSDF_DISPERSIVE))
				{
					state.chromatic=false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					pcol *= wl_col;
				}

				ray.from = sp.P;
				ray.dir = wo;
				ray.tmin = MIN_RAYDIST;
				ray.tmax = -1.0;
				++nBounces;
			}
		}

		pb->done();
		pb->setTag("Caustics photon map built.");
		delete lightPowerD;
		
		Y_INFO << integratorName << ": Shot "<<curr<<" caustic photons from " << numCLights <<" light(s)." << yendl;
		Y_VERBOSE << integratorName << ": Stored caustic photons: " << session.causticMap->nPhotons() << yendl;
	}
	else
	{
		Y_INFO << integratorName << ": No caustic source lights found, skiping caustic gathering..." << yendl;
	}

	Y_INFO << integratorName << ": Shot "<<curr<<" caustic photons from " << numCLights <<" light(s)." << yendl;
	Y_INFO << integratorName << ": Stored caustic photons: " << causticMap.nPhotons() << yendl;
	Y_INFO << integratorName << ": Stored diffuse photons: " << diffuseMap.nPhotons() << yendl;

	if(diffuseMap.nPhotons() > 0)
	{
		Y_INFO << integratorName << ": Building diffuse photons kd-tree:" << yendl;
		pb->setTag("Building diffuse photons kd-tree...");
		diffuseMap.updateTree();
		Y_INFO << integratorName << ": Done." << yendl;
	}
		
	tmplights.clear();

	std::thread * causticMapBuildKdTree_thread = nullptr;
	
	if(usePhotonCaustics && session.causticMap->nPhotons() > 0 && scene->getNumThreadsPhotons() >= 2)
	{
		Y_INFO << integratorName << ": Building caustic photons kd-tree:" << yendl;
		pb->setTag("Building caustic photons kd-tree...");

		causticMapBuildKdTree_thread = new std::thread(&photonIntegrator_t::photonMapKdTreeWorker, this, session.causticMap);
	}
	else
	{
		if( usePhotonCaustics && session.causticMap->nPhotons() > 0)
		{
			Y_INFO << integratorName << ": Building caustic photons kd-tree:" << yendl;
			pb->setTag("Building caustic photons kd-tree...");
			session.causticMap->updateTree();
			Y_VERBOSE << integratorName << ": Done." << yendl;
		}
	}
	// SSS --------->
    if (usePhotonSSS)
	{
        //Y_INFO << "SSSMap : " << SSSMaps.size() << yendl;
        //createSSSMaps();
        // prepass
        createSSSMapsByPhotonTracing();
        // povman:add text to banner info
		set << "SSS shoot:" << nSSSPhotons << " photons. ";
		std::map<const object3d_t*, photonMap_t*>::iterator it = SSSMaps.begin();
		while (it!=SSSMaps.end())
        {
            // change to use own updateTree function
		    it->second->updatePhTree();
            Y_INFO << integratorName << ": Building SSS photons kd-tree:" << yendl;
            it++;
        }
	}
	// end

	if( usePhotonDiffuse && session.diffuseMap->nPhotons() > 0 && scene->getNumThreadsPhotons() >= 2 && diffuseMapBuildKdTree_thread)
	{
		diffuseMapBuildKdTree_thread->join();
		delete diffuseMapBuildKdTree_thread;
		diffuseMapBuildKdTree_thread = nullptr;

		Y_VERBOSE << integratorName << ": Diffuse photon map: done." << yendl;
	}

	lookupRad = 4*dsRadius*dsRadius;

	tmplights.clear();

	if(!intpb) delete pb;

	if(finalGather) //create radiance map:
	{
		// == remove too close radiance points ==//
        kdtree::pointKdTree< radData_t > *rTree = new kdtree::pointKdTree< radData_t >(pgdat.rad_points);
		std::vector< radData_t > cleaned;
		for(unsigned int i=0; i<pgdat.rad_points.size(); ++i)
		{
			if(pgdat.rad_points[i].use)
			{
				cleaned.push_back(pgdat.rad_points[i]);
				eliminatePhoton_t elimProc(pgdat.rad_points[i].normal);
				float maxrad = 0.01f*dsRadius; // 10% of diffuse search radius
				rTree->lookup(pgdat.rad_points[i].pos, elimProc, maxrad);
			}
		}		
		pgdat.rad_points.swap(cleaned);
		// ================ //
		int nThreads = scene->getNumThreads();
		pgdat.radianceVec.resize(pgdat.rad_points.size());
		if(intpb) pgdat.pbar = intpb;
		else pgdat.pbar = new ConsoleProgressBar_t(80);
		pgdat.pbar->init(pgdat.rad_points.size());
		pgdat.pbar->setTag("Pregathering radiance data for final gathering...");
		std::vector<preGatherWorker_t *> workers;
		for(int i=0; i<nThreads; ++i) workers.push_back(new preGatherWorker_t(&pgdat, dsRadius, nDiffuseSearch));

		for(int i=0;i<nThreads;++i) workers[i]->run();
		for(int i=0;i<nThreads;++i)	workers[i]->wait();
		for(int i=0;i<nThreads;++i)	delete workers[i];

		radianceMap.swapVector(pgdat.radianceVec);
		pgdat.pbar->done();
		pgdat.pbar->setTag("Pregathering radiance data done...");
		if(!intpb) delete pgdat.pbar;
		Y_VERBOSE << integratorName << ": Radiance tree built... Updating the tree..." << yendl;
		session.radianceMap->updateTree();
		Y_VERBOSE << integratorName << ": Done." << yendl;
		
		delete rTree;
		rTree = nullptr;
	}

	if(usePhotonCaustics && session.causticMap->nPhotons() > 0 && scene->getNumThreadsPhotons() >= 2 && causticMapBuildKdTree_thread)
	{
		causticMapBuildKdTree_thread->join();
		delete causticMapBuildKdTree_thread;
		causticMapBuildKdTree_thread = nullptr;

		Y_VERBOSE << integratorName << ": Caustic photon map: done." << yendl;
	}

	if(photonMapProcessing == PHOTONS_GENERATE_AND_SAVE)
	{
		if( usePhotonDiffuse )
		{
			pb->setTag("Saving diffuse photon map to file...");
			std::string filename = session.getPathImageOutput() + "_diffuse.photonmap";
			Y_INFO << integratorName << ": Saving diffuse photon map to: " << filename << yendl;
			if(photonMapSave(session.diffuseMap, filename)) Y_VERBOSE << integratorName << ": Diffuse map saved." << yendl;
		}

		Y_INFO << integratorName << ": Creating radiance map..." << yendl;
		progressBar_t *pbar;
		if(intpb) pbar = intpb;
		else pbar = new ConsoleProgressBar_t(80);
		pbar->init(pgdat.rad_points.size());
		foundPhoton_t *gathered = (foundPhoton_t *)malloc(nDifuseSearch * sizeof(foundPhoton_t));
		PFLOAT dsRadius_2 = dsRadius*dsRadius;
		for(unsigned int n=0; n< pgdat.rad_points.size(); ++n)
		{
			PFLOAT radius = dsRadius_2; //actually the square radius...
			int nGathered = diffuseMap.gather(pgdat.rad_points[n].pos, gathered, nDifuseSearch, radius);
			color_t sum(0.0);
			if(nGathered > 0)
			{
				color_t surfCol = pgdat.rad_points[n].refl;
				vector3d_t rnorm = pgdat.rad_points[n].normal;
				float scale = 1.f / ( float(diffuseMap.nPaths()) * radius * M_PI);

				if(isnan(scale))
				{
					Y_WARNING << integratorName << ": NaN on (scale)" << yendl;
					break;
				}

				for(int i=0; i<nGathered; ++i)
				{
					vector3d_t pdir = gathered[i].photon->direction();

					if( rnorm * pdir > 0.f ) sum += surfCol * scale * gathered[i].photon->color();
					else sum += pgdat.rad_points[n].transm * scale * gathered[i].photon->color();
				}
			}
			photon_t radP(pgdat.rad_points[n].normal, pgdat.rad_points[n].pos, sum);
			radianceMap.pushPhoton(radP);
			if(n && !(n&7)) pbar->update(8);
		}
	}

	gTimer.stop("prepass");
	Y_INFO << integratorName << ": Photonmap building time: " << std::fixed << std::setprecision(1) << gTimer.getTime("prepass") << "s" << " (" << scene->getNumThreadsPhotons() << " thread(s))" << yendl;

	set << "| photon maps: " << std::fixed << std::setprecision(1) << gTimer.getTime("prepass") << "s" << " [" << scene->getNumThreadsPhotons() << " thread(s)]";

	yafLog.appendRenderSettings(set.str());
	
	for (std::string line; std::getline(set, line, '\n');) Y_VERBOSE << line << yendl;

	return true;
}

// final gathering: this is basically a full path tracer only that it uses the radiance map only
// at the path end. I.e. paths longer than 1 are only generated to overcome lack of local radiance detail.
// precondition: initBSDF of current spot has been called!
color_t photonIntegrator_t::finalGathering(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, colorPasses_t &colorPasses) const
{
	color_t pathCol(0.0);
	void *first_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE+7];
	void *n_udat = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	const volumeHandler_t *vol;
	color_t vcol(0.f);
	float W = 0.f;

	int nSampl = std::max(1, nPaths/state.rayDivision);
	for(int i=0; i<nSampl; ++i)
	{
		color_t throughput( 1.0 );
		float length=0;
		surfacePoint_t hit=sp;
		vector3d_t pwo = wo;
		ray_t pRay;
		BSDF_t matBSDFs;
		bool did_hit;
		const material_t *p_mat = sp.material;
		unsigned int offs = nPaths * state.pixelSample + state.samplingOffs + i; // some redundancy here...
		color_t lcol, scol;
		// "zero'th" FG bounce:
		float s1 = RI_vdC(offs);
		float s2 = scrHalton(2, offs);
		if(state.rayDivision > 1)
		{
			s1 = addMod1(s1, state.dc1);
			s2 = addMod1(s2, state.dc2);
		}

		sample_t s(s1, s2, BSDF_DIFFUSE|BSDF_REFLECT|BSDF_TRANSMIT); // glossy/dispersion/specular done via recursive raytracing
		scol = p_mat->sample(state, hit, pwo, pRay.dir, s, W);

		scol *= W;
		if(scol.isBlack()) continue;

		pRay.tmin = scene->rayMinDist;
		pRay.tmax = -1.0;
		pRay.from = hit.P;
		throughput = scol;

		if( !(did_hit = scene->intersect(pRay, hit)) ) continue; //hit background

		p_mat = hit.material;
		length = pRay.tmax;
		state.userdata = n_udat;
		matBSDFs = p_mat->getFlags();
		bool has_spec = matBSDFs & BSDF_SPECULAR;
		bool caustic = false;
		bool close = length < gatherDist;
		bool do_bounce = close || has_spec;
		// further bounces construct a path just as with path tracing:
		for(int depth=0; depth<gatherBounces && do_bounce; ++depth)
		{
			int d4 = 4*depth;
			pwo = -pRay.dir;
			p_mat->initBSDF(state, hit, matBSDFs);

			if((matBSDFs & BSDF_VOLUMETRIC) && (vol=p_mat->getVolumeHandler(hit.N * pwo < 0)))
			{
				if(vol->transmittance(state, pRay, vcol)) throughput *= vcol;
			}

			if(matBSDFs & (BSDF_DIFFUSE))
			{
				if(close)
				{
					lcol = estimateOneDirectLight(state, hit, pwo, offs, tmpColorPasses);
				}
				else if(caustic)
				{
					vector3d_t sf = FACE_FORWARD(hit.Ng, hit.N, pwo);
					const photon_t *nearest = session.radianceMap->findNearest(hit.P, sf, lookupRad);
					if(nearest) lcol = nearest->color();
				}

				if(close || caustic)
				{
					if(matBSDFs & BSDF_EMIT) lcol += p_mat->emit(state, hit, pwo);
					pathCol += lcol*throughput;
				}
			}

			s1 = scrHalton(d4+3, offs);
			s2 = scrHalton(d4+4, offs);

			if(state.rayDivision > 1)
			{
				s1 = addMod1(s1, state.dc1);
				s2 = addMod1(s2, state.dc2);
			}

			sample_t sb(s1, s2, (close) ? BSDF_ALL : BSDF_ALL_SPECULAR | BSDF_FILTER);
			scol = p_mat->sample(state, hit, pwo, pRay.dir, sb, W);

			if( sb.pdf <= 1.0e-6f)
			{
				did_hit=false;
				break;
			}

			scol *= W;

			pRay.tmin = scene->rayMinDist;
			pRay.tmax = -1.0;
			pRay.from = hit.P;
			throughput *= scol;
			did_hit = scene->intersect(pRay, hit);

			if(!did_hit) //hit background
			{
				 if(caustic && background && background->hasIBL() && background->shootsCaustic())
				 {
					pathCol += throughput * (*background)(pRay, state, true);
				 }
				 break;
			}

			p_mat = hit.material;
			length += pRay.tmax;
			caustic = (caustic || !depth) && (sb.sampledFlags & (BSDF_SPECULAR | BSDF_FILTER));
			close =  length < gatherDist;
			do_bounce = caustic || close;
		}

		if(did_hit)
		{
			p_mat->initBSDF(state, hit, matBSDFs);
			if(matBSDFs & (BSDF_DIFFUSE | BSDF_GLOSSY))
			{
				vector3d_t sf = FACE_FORWARD(hit.Ng, hit.N, -pRay.dir);
				const photon_t *nearest = session.radianceMap->findNearest(hit.P, sf, lookupRad);
				if(nearest) lcol = nearest->color();
				if(matBSDFs & BSDF_EMIT) lcol += p_mat->emit(state, hit, -pRay.dir);
				pathCol += lcol * throughput;
			}
		}
		state.userdata = first_udat;
	}
	return pathCol / (float)nSampl;
}

colorA_t photonIntegrator_t::integrate(renderState_t &state, diffRay_t &ray, colorPasses_t &colorPasses, int additionalDepth /*=0*/) const
{
	static int _nMax=0;
	static int calls=0;
	++calls;
	color_t col(0.0);
	CFLOAT alpha = 1.0;
	surfacePoint_t sp;

	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;

	if (transpBackground) alpha = 0.0;

	if(scene->intersect(ray, sp))
	{
		unsigned char userdata[USER_DATA_SIZE+7];
		state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
		if(state.raylevel == 0)
		{
			state.chromatic = true;
			state.includeLights = true;
		}
		BSDF_t bsdfs;
		//vector3d_t N_nobump = sp.N; // unused
		vector3d_t wo = -ray.dir;
		const material_t *material = sp.material;
		material->initBSDF(state, sp, bsdfs);

		if(additionalDepth < material->getAdditionalDepth()) additionalDepth = material->getAdditionalDepth();
		
		col += colorPasses.probe_add(PASS_INT_EMIT, material->emit(state, sp, wo), state.raylevel == 0);
		
		state.includeLights = false;
		spDifferentials_t spDiff(sp, ray);

		if(finalGather)
		{
			if(showMap)
			{
				vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
				const photon_t *nearest = session.radianceMap->findNearest(sp.P, N, lookupRad);
				if(nearest) col += nearest->color();
			}
			else
			{
				if(state.raylevel == 0 && colorPasses.enabled(PASS_INT_RADIANCE))
				{
					vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
					const photon_t *nearest = session.radianceMap->findNearest(sp.P, N, lookupRad);
					if(nearest) colorPasses(PASS_INT_RADIANCE) = nearest->color();
				}
				
				// contribution of light emitting surfaces
				if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);

				if(bsdfs & BSDF_DIFFUSE)
				{
					col += estimateAllDirectLight(state, sp, wo, colorPasses);;
					
					if(AA_clamp_indirect>0.f)
					{
						color_t tmpCol = finalGathering(state, sp, wo, colorPasses);
						tmpCol.clampProportionalRGB(AA_clamp_indirect);
						col += colorPasses.probe_set(PASS_INT_DIFFUSE_INDIRECT, tmpCol, state.raylevel == 0);
					}
					else col += colorPasses.probe_set(PASS_INT_DIFFUSE_INDIRECT, finalGathering(state, sp, wo, colorPasses), state.raylevel == 0);
				}
			}
		}
		else
		{
			if(usePhotonDiffuse && showMap)
			{
				vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
				const photon_t *nearest = session.diffuseMap->findNearest(sp.P, N, dsRadius);
				if(nearest) col += nearest->color();
			}
			else
			{
				if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);

				if(bsdfs & BSDF_DIFFUSE)
				{
					col += estimateAllDirectLight(state, sp, wo, colorPasses);
				}
				
				foundPhoton_t *gathered = (foundPhoton_t *)alloca(nDiffuseSearch * sizeof(foundPhoton_t));
				float radius = dsRadius; //actually the square radius...

				int nGathered=0;

				if(diffuseMap.nPhotons() > 0) nGathered = diffuseMap.gather(sp.P, gathered, nDiffuseSearch, radius);
				color_t sum(0.0);
				if(usePhotonDiffuse && nGathered > 0)
				{
					if(nGathered > _nMax) _nMax = nGathered;

					float scale = 1.f / ( (float)session.diffuseMap->nPaths() * radius * M_PI);
					for(int i=0; i<nGathered; ++i)
					{
						vector3d_t pdir = gathered[i].photon->direction();
						color_t surfCol = material->eval(state, sp, wo, pdir, BSDF_DIFFUSE);

						col += colorPasses.probe_add(PASS_INT_DIFFUSE_INDIRECT, surfCol * scale * gathered[i].photon->color(), state.raylevel == 0);
					}
				}
			}
		}
		// add caustics
		if(bsdfs & BSDF_DIFFUSE) col += estimateCausticPhotons(state, sp, wo);

        // povman: if have translucent SSS material ..
        if(bsdfs & BSDF_TRANSLUCENT)
		{
            // ..and 'use SSS photons' is active.
            if(usePhotonSSS)
            {
                col += estimateSSSMaps(state, sp, wo); //(1)
                col += estimateSSSSingleSImportantSampling(state, sp, wo);
            }
		}
		//end
		recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);

			if(colorPasses.enabled(PASS_INT_AO_CLAY))
			{
				colorPasses(PASS_INT_AO_CLAY) = sampleAmbientOcclusionPassClay(state, sp, wo);
			}
		}
		
		if(transpRefractedBackground)
		{
			float m_alpha = material->getAlpha(state, sp, wo);
			alpha = m_alpha + (1.f-m_alpha)*alpha;
		}
		else alpha = 1.0;
	}
	else //nothing hit, return background
	{
		if(background && !transpRefractedBackground)
		{
			col += colorPasses.probe_set(PASS_INT_ENV, (*background)(ray, state), state.raylevel == 0);
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;
    //
    colorA_t result = col;
    colorA_t colVolTransmitt = scene->volIntegrator->transmittance(state, ray);
    result *= colVolTransmitt;
    result += scene->volIntegrator->integrate(state, ray);
	result.A = std::max(alpha, 1.f - colVolTransmitt.A);

    return result;
}

integrator_t* photonIntegrator_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false;
	bool finalGather=true;
	bool show_map=false;
	int shadowDepth=5;
	int raydepth=5;
	int numPhotons = 100000;
	int numCPhotons = 500000;
	int search = 50;
	int caustic_mix = 50;
	int bounces = 5;
	int fgPaths = 32;
	int fgBounces = 2;
	float dsRad=0.1;
	float cRad=0.01;
	float gatherDist=0.2;
	bool bg_transp = true;
	bool bg_transp_refract = true;

	// SSS
    bool useSSS=false;
	int sssdepth = 10, sssPhotons = 200000;
	int singleSSamples = 128;
	float sScale = 40.f;


	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("raydepth", raydepth);
	params.getParam("photons", numPhotons);
	params.getParam("cPhotons", numCPhotons);
	params.getParam("diffuseRadius", dsRad);
	params.getParam("causticRadius", cRad);
	params.getParam("search", search);
	caustic_mix = search;
	params.getParam("caustic_mix", caustic_mix);
	params.getParam("bounces", bounces);
	params.getParam("finalGather", finalGather);
	params.getParam("fg_samples", fgPaths);
	params.getParam("fg_bounces", fgBounces);
	gatherDist = dsRad;
	params.getParam("fg_min_pathlen", gatherDist);
	params.getParam("show_map", show_map);
	params.getParam("bg_transp", bg_transp);
	params.getParam("bg_transp_refract", bg_transp_refract);
	// SSS
	params.getParam("useSSS", useSSS);
	params.getParam("sssPhotons", sssPhotons);
	params.getParam("sssDepth", sssdepth);
	params.getParam("singleScatterSamples", singleSSamples);
	params.getParam("sssScale", sScale);

	photonIntegrator_t* ite = new photonIntegrator_t(numPhotons, numCPhotons, transpShad, shadowDepth, dsRad, cRad);
	
	ite->usePhotonCaustics = caustics;
	ite->usePhotonDiffuse = diffuse;
	
	ite->rDepth = raydepth;
	ite->nDiffuseSearch = search;
	ite->nCausSearch = caustic_mix;
	ite->finalGather = finalGather;
	ite->maxBounces = bounces;
	ite->causDepth = bounces;
	ite->nPaths = fgPaths;
	ite->gatherBounces = fgBounces;
	ite->showMap = show_map;
	ite->gatherDist = gatherDist;
	// Background settings
	ite->transpBackground = bg_transp;
	ite->transpRefractedBackground = bg_transp_refract;
	// SSS
	ite->usePhotonSSS = useSSS;
	ite->nSSSPhotons = sssPhotons;
	ite->nSSSDepth = sssdepth;
	ite->nSingleScatterSamples = singleSSamples;
	ite->isDirectLight = false;
	ite->sssScale = sScale;

	return ite;
}

extern "C"
{
	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("photonmapping", photonIntegrator_t::factory);
	}
}

__END_YAFRAY
