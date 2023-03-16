/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifdef WIN32
#include <windows.h>
#include <Psapi.h>
#endif

#include <QDir>
#include <QTemporaryDir>
#include <QTemporaryFile>

#include <thread>
#include <vector>

#include "filter_iPSR.h"
#include "poisson_utils.h"
#include "Src/kdtree.h"

using namespace std;

FilteriPSRPlugin::FilteriPSRPlugin()
{
	typeList = {FP_SCREENED_POISSON};

	for (ActionIDType tt : types()){
		actionList.push_back(new QAction(filterName(tt), this));
	}
}

FilteriPSRPlugin::~FilteriPSRPlugin()
{
}

QString FilteriPSRPlugin::pluginName() const
{
	return "FilteriPSR";
}

QString FilteriPSRPlugin::filterName(ActionIDType filter) const
{
	if (filter == FP_SCREENED_POISSON) {
		return "Iterative Poisson Surface Reconstruction";
	}
	else {
		assert(0);
		return QString();
	}
}

QString FilteriPSRPlugin::pythonFilterName(ActionIDType f) const
{
	if (f == FP_SCREENED_POISSON) {
		return "iterative_Poisson_Surface_Reconstruction";
	}
	else {
		assert(0);
		return QString();
	}
}

QString FilteriPSRPlugin::filterInfo(ActionIDType filter) const
{
	if (filter == FP_SCREENED_POISSON)
		return	"iPSR extends the popular Poisson Surface Reconstruction (https://github.com/mkazhdan/PoissonRecon). iPSR has no more need of oriented normals as input, but infers the normals in an iterative manner. It is used to reconstruct surface from only points input.";
	else {
		return "Error!";
	}
}

FilterPlugin::FilterClass FilteriPSRPlugin::getClass(const QAction* a) const
{
	if (ID(a) == FP_SCREENED_POISSON){
		return FilterPlugin::Remeshing;
	}
	else {
		assert(0);
		return FilterPlugin::Generic;
	}
}

int FilteriPSRPlugin::getRequirements(const QAction* a)
{
	if (ID(a) == FP_SCREENED_POISSON) {
		return MeshModel::MM_NONE;
	}
	else {
		assert(0);
		return 0;
	}
}

std::map<std::string, QVariant> FilteriPSRPlugin::applyFilter(
		const QAction* filter,
		const RichParameterList& params,
		MeshDocument& md,
		unsigned int& /*postConditionMask*/,
		vcg::CallBackPos* cb)
{
	bool currDirChanged=false;
	QDir currDir = QDir::current();

	if (ID(filter) == FP_SCREENED_POISSON) {
		//Using tmp dir
		QTemporaryDir tmpdir;
		QTemporaryFile file(tmpdir.path());
		if (!file.open()) { //if a file cannot be created in the tmp folder
			log("Warning - tmp folder is not writable.");

			QTemporaryFile file2("./_tmp_XXXXXX.tmp"); //try to create a file in the meshlab folder
			if (!file2.open()){ //if a file cannot be created in the tmp and in the meshlab folder, we cannot run the filter
				log("Warning - current folder is not writable. Screened Poisson Merging needs to save intermediate files in the tmp working folder. Project and meshes must be in a write-enabled folder. Please save your data in a suitable folder before applying.");
				throw MLException("current and tmp folder are not writable.<br> Screened Poisson Merging needs to save intermediate files in the current working folder.<br> Project and meshes must be in a write-enabled folder.<br> Please save your data in a suitable folder before applying.");
			}
		}
		else { //if the tmp folder is writable, we will use it
			currDirChanged=true;
			QDir::setCurrent(tmpdir.path());
		}

		PoissonParam<Scalarm> pp;
		pp.MaxDepthVal = params.getInt("depth");
		pp.FullDepthVal = params.getInt("fullDepth");
		pp.CGDepthVal= params.getInt("cgDepth");
		pp.ScaleVal = params.getFloat("scale");
		pp.SamplesPerNodeVal = params.getFloat("samplesPerNode");
		pp.PointWeightVal = params.getFloat("pointWeight");
		pp.ItersVal = params.getInt("iters");
		pp.ConfidenceFlag = params.getBool("confidence");
		pp.DensityFlag = true;
		pp.CleanFlag = params.getBool("preClean");
		pp.ThreadsVal = params.getInt("threads");

		bool goodNormal=true, goodColor=true;
		if(params.getBool("visibleLayer") == false) {
			PoissonClean(md.mm()->cm, pp.ConfidenceFlag, pp.CleanFlag);
			goodNormal=HasGoodNormal(md.mm()->cm);
			goodColor = md.mm()->hasDataMask(MeshModel::MM_VERTCOLOR);
		}
		else {
			MeshModel *_mm=md.nextVisibleMesh();
			while(_mm != nullptr) {
				PoissonClean(_mm->cm,  pp.ConfidenceFlag, pp.CleanFlag);
				goodNormal &= HasGoodNormal(_mm->cm);
				goodColor  &= _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
				_mm=md.nextVisibleMesh(_mm);
			}
		}
/*
		if(!goodNormal) {
			throw MLException("Filter requires correct per vertex normals.<br>"
								 "E.g. it is necessary that your <b>ALL</b> the input vertices have a proper, not-null normal.<br> "
								 "Try enabling the <i>pre-clean<i> option and retry.<br><br>"
								 "To permanently remove this problem:<br>"
								 "If you encounter this error on a triangulated mesh try to use the <i>Remove Unreferenced Vertices</i> filter"
								 "If you encounter this error on a pointcloud try to use the <i>Conditional Vertex Selection</i> filter"
								 "with function '(nx==0.0) && (ny==0.0) && (nz==0.0)', and then <i>delete selected vertices</i>.<br>");
		}
*/	
		//randomise normals
		typedef double Real;
		Point3D<Real> zero_normal(Point3D<Real>(0,0,0));
		srand(0);
		CMeshO &m = md.mm()->cm;
		for (size_t i = 0; i <md.mm()->cm.vn; i++) {
			m.vert[i].N()[0]=rand() % 1001 - 500.0, rand() % 1001 - 500.0;
			m.vert[i].N()[1]=rand() % 1001 - 500.0, rand() % 1001 - 500.0;
			m.vert[i].N()[2]=rand() % 1001 - 500.0, rand() % 1001 - 500.0;		
			//still need to normalise
			Point3D<Real> n(m.vert[i].N()[0], m.vert[i].N()[1], m.vert[i].N()[2]);
			while(n[0] == zero_normal[0] && n[1] == zero_normal[1] && n[2] == zero_normal[2]);
				m.vert[i].N().Normalize();
		}
		

		

		// construct the Kd-Tree
		kdt::KDTree<kdt::KDTreePoint> tree;
		{vector<kdt::KDTreePoint> vertices;
			vertices.reserve(m.vn);
			for (size_t i = 0; i < m.vn; i++) {
				array<double,3> p = {m.vert[i].P()[0], m.vert[i].P()[1], m.vert[i].P()[2]};
				vertices.push_back(kdt::KDTreePoint(p));
			}
			tree.build(vertices);
		}
		int epoch =0;
		int iters =params.getInt("iters");
		float pointweight =params.getFloat("pointWeight");
		int depth =params.getInt("depth");
		int k_neighbors =params.getInt("neighbors");

		while (epoch<iters)
		{	
			++epoch;
			// Poisson reconstruction
			log("Iter: %d\n", epoch);

			MeshModel *pm =md.addNewMesh("","temp iPSR Mesh",false);
			
			
			MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
			_Execute<Scalarm,2,BOUNDARY_NEUMANN,PlyColorAndValueVertex<Scalarm> >(&meshStream,md.mm()->cm.bbox,pm->cm,pp,cb);
			CMeshO &mn = md.getMesh(pm -> id())->cm;

			vector<vector<int>> nearestSamples(mn.fn);
			vector<Point3D<Real>> normals(mn.fn);
			log("Number of faces: %d\n", mn.fn);
			for (int i = 0 ; i < mn.fn ;i++)
			{	
				Point3D<Real> c0;
				Point3D<Real> c1;
				Point3D<Real> c2;
				Point3D<Real> c;
				//get Point of each face
				for (int j = 0; j < 3; j++) {
					c0[j] = mn.face[i].V(0) -> P()[j];
					c1[j] = mn.face[i].V(1) -> P()[j];
					c2[j] = mn.face[i].V(2) -> P()[j];
				}
				c = c0 + c1 + c2;

				c /= 3.0;
				array<Real,3> a {c[0], c[1], c[2]};
				nearestSamples[i] = tree.knnSearch(kdt::KDTreePoint(a), k_neighbors);
				Point3D<Real> n;
				CrossProduct(c1-c0, c2-c0,n);
				normals[i] =n;
				
			}
			// update sample point normals
			vector<Point3D<Real>> projective_normals(mn.vn, zero_normal);
			for (size_t i = 0; i < nearestSamples.size(); i++)
			{
				for (size_t j = 0; j < nearestSamples[i].size(); ++j)
				{
				projective_normals[nearestSamples[i][j]][0] += normals[i][0];
				projective_normals[nearestSamples[i][j]][1] += normals[i][1];
				projective_normals[nearestSamples[i][j]][2] += normals[i][2];
				}
			}
			//still need to normilize
			//for (int i = 0; i < (int)projective_normals.size(); ++i)
			//	projective_normals[i].normalize();
			
			// compute the average normal variation of the top 1/1000 points
			size_t heap_size = static_cast<size_t>(ceil(mn.vn / 1000.0));
			priority_queue<double, vector<double>, greater<double>> min_heap;
			for (size_t i = 0; i < m.vn; ++i)
			{
				if (!(projective_normals[i][0] == zero_normal[0]&& projective_normals[i][1] == zero_normal[1] && projective_normals[i][2] == zero_normal[2]))
				{
					Point3D<Real> normal;
					normal[0] = m.vert[i].N()[0];
					normal[1] = m.vert[i].N()[1];
					normal[2] = m.vert[i].N()[2];
					double diff = Point3D<Real>::Dot(projective_normals[i]- normal,projective_normals[i]- normal);
					if (min_heap.size() < heap_size)
						min_heap.push(diff);
					else if (diff > min_heap.top())
					{
						min_heap.pop();
						min_heap.push(diff);
					}
					for (size_t j = 0; j < 3; ++j)
						m.vert[i].N()[j]=projective_normals[i][j];
				}
			}

			heap_size = min_heap.size();
			double ave_max_diff = 0;
			while (!min_heap.empty())
			{
				ave_max_diff += sqrt(min_heap.top());
				min_heap.pop();
			}
			ave_max_diff /= heap_size;
			log("normals variation %f\n", ave_max_diff);
			md.delMesh(pm->id());
			if (ave_max_diff < 0.100075)
				break;



		}
		
		MeshModel *pm =md.addNewMesh("","iPSR Mesh",false);

		//md.setVisible(pm->id(),true);
		pm->updateDataMask(MeshModel::MM_VERTQUALITY);
		
		if(goodColor)
			pm->updateDataMask(MeshModel::MM_VERTCOLOR);

		MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
		_Execute<Scalarm,2,BOUNDARY_NEUMANN,PlyColorAndValueVertex<Scalarm> >(&meshStream,md.mm()->cm.bbox,pm->cm,pp,cb);
		
		pm->updateBoxAndNormals();
		md.setVisible(pm->id(),true);
		md.setCurrentMesh(pm->id());
		if(currDirChanged)
			QDir::setCurrent(currDir.path());
	}
	else {
		wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}

RichParameterList FilteriPSRPlugin::initParameterList(
		const QAction* filter,
		const MeshModel&)
{
	RichParameterList parlist;
	unsigned int nThreads = std::thread::hardware_concurrency();
	if (nThreads == 0) nThreads = 8;
	if (ID(filter) == FP_SCREENED_POISSON) {
		parlist.addParam(RichInt ("iters",
                                            30,
                                            "maximum number of iterations",
                                            "The maximum number of iterations. The default value of this parameter is 30.\n"));
            parlist.addParam(RichFloat ("pointWeight",
                                            10,
                                            "interpolation weight",
                                            "The pointWeight parameter of screened Poisson surface reconstruction. The default value for this parameter is 10.\n"));
            parlist.addParam(RichInt ("depth",
                                            10,
                                            "Samples per Node",
                                            "The depth parameter of screened Poisson surface reconstruction. It is the maximum possible depth of the octree. The default value of this parameter is 10.\n"));
            parlist.addParam(RichInt ("neighbors",
                                             10,
                                             "number of neighbors",
                                             "The number of the closest sample points to search from every reconstructed triangle face. The suggested range is between 10 and 20. The default value of this parameter is 10.\n"));

		parlist.addParam(RichBool("visibleLayer", false, "Merge all visible layers", "Enabling this flag means that all the visible layers will be used for providing the points."));
		//parlist.addParam(RichInt("depth", 8, "Reconstruction Depth", "This integer is the maximum depth of the tree that will be used for surface reconstruction. Running at depth d corresponds to solving on a voxel grid whose resolution is no larger than 2^d x 2^d x 2^d. Note that since the reconstructor adapts the octree to the sampling density, the specified reconstruction depth is only an upper bound. The default value for this parameter is 8."));
		parlist.addParam(RichInt("fullDepth", 5, "Adaptive Octree Depth", "This integer specifies the depth beyond depth the octree will be adapted. At coarser depths, the octree will be complete, containing all 2^d x 2^d x 2^d nodes. The default value for this parameter is 5.", true));
		parlist.addParam(RichInt("cgDepth", 0, "Conjugate Gradients Depth", "This integer is the depth up to which a conjugate-gradients solver will be used to solve the linear system. Beyond this depth Gauss-Seidel relaxation will be used. The default value for this parameter is 0.", true));
		parlist.addParam(RichFloat("scale", 1.1, "Scale Factor", "This floating point value specifies the ratio between the diameter of the cube used for reconstruction and the diameter of the samples' bounding cube. The default value is 1.1.", true));
		parlist.addParam(RichFloat("samplesPerNode", 1.5, "Minimum Number of Samples", "This floating point value specifies the minimum number of sample points that should fall within an octree node as the octree construction is adapted to sampling density. For noise-free samples, small values in the range [1.0 - 5.0] can be used. For more noisy samples, larger values in the range [15.0 - 20.0] may be needed to provide a smoother, noise-reduced, reconstruction. The default value is 1.5."));
		//parlist.addParam(RichFloat("pointWeight", 4, "Interpolation Weight", "This floating point value specifies the importants that interpolation of the point samples is given in the formulation of the screened Poisson equation. The results of the original (unscreened) Poisson Reconstruction can be obtained by setting this value to 0. The default value for this parameter is 4."));
		//parlist.addParam(RichInt("iters", 8, "Gauss-Seidel Relaxations", "This integer value specifies the number of Gauss-Seidel relaxations to be performed at each level of the hierarchy. The default value for this parameter is 8.", true));
		parlist.addParam(RichBool("confidence", false, "Confidence Flag", "Enabling this flag tells the reconstructor to use the quality as confidence information; this is done by scaling the unit normals with the quality values. When the flag is not enabled, all normals are normalized to have unit-length prior to reconstruction."));
		parlist.addParam(RichBool("preClean", false, "Pre-Clean", "Enabling this flag force a cleaning pre-pass on the data removing all unreferenced vertices or vertices with null normals."));
		parlist.addParam(RichInt("threads", nThreads, "Number Threads", "Maximum number of threads that the reconstruction algorithm can use."));
	}
	return parlist;
}

int FilteriPSRPlugin::postCondition(const QAction* filter) const
{
	if (ID(filter) == FP_SCREENED_POISSON){
		return MeshModel::MM_VERTNUMBER + MeshModel::MM_FACENUMBER;
	}
	else {
		return MeshModel::MM_ALL;
	}
}


FilterPlugin::FilterArity FilteriPSRPlugin::filterArity(const QAction*) const
{
	return VARIABLE;
}

