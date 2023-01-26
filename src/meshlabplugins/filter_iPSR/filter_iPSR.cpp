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

//#include <QDir>
//#include <QTemporaryDir>
//#include <QTemporaryFile>

#include <thread>



#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include "src/kdtree.h"
#include "src/utility.h"
#include "src/PoissonRecon.h"
#include "src/PointStreamData.h"




#include "filter_iPSR.h"
//#include "poisson_utils.h"


using namespace std;
FilterScreenedPoissonPlugin::FilterScreenedPoissonPlugin()
{
	typeList = {FP_SCREENED_POISSON};

	for (ActionIDType tt : types()){
		actionList.push_back(new QAction(filterName(tt), this));
	}
}

FilterScreenedPoissonPlugin::~FilterScreenedPoissonPlugin()
{
}

QString FilterScreenedPoissonPlugin::pluginName() const
{
	return "FilteriPSR";
}

QString FilterScreenedPoissonPlugin::filterName(ActionIDType filter) const
{
	if (filter == FP_SCREENED_POISSON) {
		return "Iterative Poisson Surface Reconstruction";
	}
	else {
		assert(0);
		return QString();
	}
}

QString FilterScreenedPoissonPlugin::pythonFilterName(ActionIDType f) const
{
	if (f == FP_SCREENED_POISSON) {
		return "iterative_Poisson_Surface_Reconstruction";
	}
	else {
		assert(0);
		return QString();
	}
}

QString FilterScreenedPoissonPlugin::filterInfo(ActionIDType filter) const
{
	if (filter == FP_SCREENED_POISSON)
		return	"iPSR extends the popular Poisson Surface Reconstruction (https://github.com/mkazhdan/PoissonRecon). iPSR has no more need of oriented normals as input, but infers the normals in an iterative manner. It is used to reconstruct surface from only points input.";
	else {
		return "Error!";
	}
}

FilterPlugin::FilterClass FilterScreenedPoissonPlugin::getClass(const QAction* a) const
{
	if (ID(a) == FP_SCREENED_POISSON){
		return FilterPlugin::Remeshing;
	}
	else {
		assert(0);
		return FilterPlugin::Generic;
	}
}

int FilterScreenedPoissonPlugin::getRequirements(const QAction* a)
{
	if (ID(a) == FP_SCREENED_POISSON) {
		return MeshModel::MM_NONE;
	}
	else {
		assert(0);
		return 0;
	}
}

std::map<std::string, QVariant> FilterScreenedPoissonPlugin::applyFilter(
		const QAction* filter,
		const RichParameterList& params,
		MeshDocument& md,
		unsigned int& /*postConditionMask*/,
		vcg::CallBackPos* cb)
{
	//parse the parameters
	int iters =params.getInt("iters");
	int pointweight =params.getInt("pointWeight");
	int depth =params.getInt("depth");
	int k_neighbors =params.getInt("neighbors");

	log("iters: %d",iters);
	log("pointweight: %d",pointweight);
	log("depth: %d",depth);
	log("k_neighbors: %d",k_neighbors);



	typedef double REAL;
	const unsigned int DIM = 3U;
    vector<pair<Point<double, 3>, NormaliSPR<double, 3>>> points_normals;
	//cancel read in, use mesh from meshlab
	//ply_reader<REAL, DIM>(input_name, points_normals);

	//fill points_normals with meshlab mesh
	NormaliSPR<REAL, DIM> n(Point<REAL, DIM>(1, 0, 0));


	for (int i = 0; i < md.mm()->cm.vert.size(); i++) {
		Point<REAL, DIM> p;
		p[0] = md.mm()->cm.vert[i].P()[0];
		p[1] = md.mm()->cm.vert[i].P()[1];
		p[2] = md.mm()->cm.vert[i].P()[2];
		points_normals.push_back(make_pair(p, n));
	}

	



	string command = "PoissonRecon --in i.ply --out o.ply --bType 2 --depth " + to_string(depth) + " --pointWeight " + to_string(pointweight);
	log("command: %s",command.c_str());
	vector<string> cmd = split(command);
	vector<char *> argv_str(cmd.size());
	for (size_t i = 0; i < cmd.size(); ++i)
		argv_str[i] = &cmd[i][0];

	XForm<REAL, DIM + 1> iXForm;
	vector<double> weight_samples;
	// sample points by the octree
	points_normals = sample_points<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, iXForm, &weight_samples);

	// initialize normals randomly
	log("random initialization...\n");
	NormaliSPR<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
	srand(0);
	for (size_t i = 0; i < points_normals.size(); ++i)
	{
		do
		{
			points_normals[i].second = Point<REAL, DIM>(rand() % 1001 - 500.0, rand() % 1001 - 500.0, rand() % 1001 - 500.0);
		} while (points_normals[i].second == zero_normal);
		normalize<REAL, DIM>(points_normals[i].second);
	}

	// construct the Kd-Tree
	kdt::KDTree<kdt::KDTreePoint> tree;
	{
		vector<kdt::KDTreePoint> vertices;
		vertices.reserve(points_normals.size());
		for (size_t i = 0; i < points_normals.size(); ++i)
		{
			array<double, 3> p{points_normals[i].first[0], points_normals[i].first[1], points_normals[i].first[2]};
			vertices.push_back(kdt::KDTreePoint(p));
		}
		tree.build(vertices);
	}

	pair<vector<Point<REAL, DIM>>, vector<vector<int>>> mesh;

	// iterations
	int epoch = 0;
	while (epoch < iters)
	{
		++epoch;
		log("Iter: %d\n", epoch);

		vector<Point<REAL, DIM>>().swap(mesh.first);
		vector<vector<int>>().swap(mesh.second);

		// Poisson reconstruction
		mesh = poisson_reconstruction<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, &weight_samples);

		vector<vector<int>> nearestSamples(mesh.second.size());
		vector<Point<REAL, DIM>> normals(mesh.second.size());

		// compute face normals and map them to sample points
#pragma omp parallel for
		for (int i = 0; i < (int)nearestSamples.size(); i++)
		{
			if (mesh.second[i].size() == 3)
			{
				Point<REAL, DIM> c = mesh.first[mesh.second[i][0]] + mesh.first[mesh.second[i][1]] + mesh.first[mesh.second[i][2]];
				c /= 3;
				array<REAL, DIM> a{c[0], c[1], c[2]};
				nearestSamples[i] = tree.knnSearch(kdt::KDTreePoint(a), k_neighbors);
				normals[i] = Point<REAL, DIM>::CrossProduct(mesh.first[mesh.second[i][1]] - mesh.first[mesh.second[i][0]], mesh.first[mesh.second[i][2]] - mesh.first[mesh.second[i][0]]);
			}
		}

		// update sample point normals
		vector<NormaliSPR<REAL, DIM>> projective_normals(points_normals.size(), zero_normal);
		for (size_t i = 0; i < nearestSamples.size(); i++)
		{
			for (size_t j = 0; j < nearestSamples[i].size(); ++j)
			{
				projective_normals[nearestSamples[i][j]].normal[0] += normals[i][0];
				projective_normals[nearestSamples[i][j]].normal[1] += normals[i][1];
				projective_normals[nearestSamples[i][j]].normal[2] += normals[i][2];
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (int)projective_normals.size(); ++i)
			normalize<REAL, DIM>(projective_normals[i]);

		// compute the average normal variation of the top 1/1000 points
		size_t heap_size = static_cast<size_t>(ceil(points_normals.size() / 1000.0));
		priority_queue<double, vector<double>, greater<double>> min_heap;
		for (size_t i = 0; i < points_normals.size(); ++i)
		{
			if (!(projective_normals[i] == zero_normal))
			{
				double diff = Point<REAL, DIM>::SquareNorm((projective_normals[i] - points_normals[i].second).normal);
				if (min_heap.size() < heap_size)
					min_heap.push(diff);
				else if (diff > min_heap.top())
				{
					min_heap.pop();
					min_heap.push(diff);
				}

				points_normals[i].second = projective_normals[i];
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
		if (ave_max_diff < 0.175)
			break;
	}

	mesh = poisson_reconstruction<REAL, DIM>((int)argv_str.size(), argv_str.data(), points_normals, &weight_samples);


	//convert mesh to meshlab format
	MeshModel *pm =md.addNewMesh("","iPSR Mesch",false);
	
	//add vertices
	for (int i = 0; i < mesh.first.size(); i++) {
		Point<REAL, DIM> p = mesh.first[i];
		//convert the point to meshlab format
		CVertexO v;
		v.P()[0] = p[0];
		v.P()[1] = p[1];
		v.P()[2] = p[2];
		//add the vertex to the mesh
		pm->cm.vert.push_back(v);
		//correct the number of vertices
		pm->cm.vn++;		
	}

	





	// bool currDirChanged=false;
	// QDir currDir = QDir::current();

	// if (ID(filter) == FP_SCREENED_POISSON) {
	// 	//Using tmp dir
	// 	QTemporaryDir tmpdir;
	// 	QTemporaryFile file(tmpdir.path());
	// 	if (!file.open()) { //if a file cannot be created in the tmp folder
	// 		log("Warning - tmp folder is not writable.");

	// 		QTemporaryFile file2("./_tmp_XXXXXX.tmp"); //try to create a file in the meshlab folder
	// 		if (!file2.open()){ //if a file cannot be created in the tmp and in the meshlab folder, we cannot run the filter
	// 			log("Warning - current folder is not writable. Screened Poisson Merging needs to save intermediate files in the tmp working folder. Project and meshes must be in a write-enabled folder. Please save your data in a suitable folder before applying.");
	// 			throw MLException("current and tmp folder are not writable.<br> Screened Poisson Merging needs to save intermediate files in the current working folder.<br> Project and meshes must be in a write-enabled folder.<br> Please save your data in a suitable folder before applying.");
	// 		}
	// 	}
	// 	else { //if the tmp folder is writable, we will use it
	// 		currDirChanged=true;
	// 		QDir::setCurrent(tmpdir.path());
	// 	}

	// 	PoissonParam<Scalarm> pp;
	// 	pp.MaxDepthVal = params.getInt("depth");
	// 	pp.FullDepthVal = params.getInt("fullDepth");
	// 	pp.CGDepthVal= params.getInt("cgDepth");
	// 	pp.ScaleVal = params.getFloat("scale");
	// 	pp.SamplesPerNodeVal = params.getFloat("samplesPerNode");
	// 	pp.PointWeightVal = params.getFloat("pointWeight");
	// 	pp.ItersVal = params.getInt("iters");
	// 	pp.ConfidenceFlag = params.getBool("confidence");
	// 	pp.DensityFlag = true;
	// 	pp.CleanFlag = params.getBool("preClean");
	// 	pp.ThreadsVal = params.getInt("threads");

	// 	bool goodNormal=true, goodColor=true;
	// 	if(params.getBool("visibleLayer") == false) {
	// 		PoissonClean(md.mm()->cm, pp.ConfidenceFlag, pp.CleanFlag);
	// 		goodNormal=HasGoodNormal(md.mm()->cm);
	// 		goodColor = md.mm()->hasDataMask(MeshModel::MM_VERTCOLOR);
	// 	}
	// 	else {
	// 		MeshModel *_mm=md.nextVisibleMesh();
	// 		while(_mm != nullptr) {
	// 			PoissonClean(_mm->cm,  pp.ConfidenceFlag, pp.CleanFlag);
	// 			goodNormal &= HasGoodNormal(_mm->cm);
	// 			goodColor  &= _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
	// 			_mm=md.nextVisibleMesh(_mm);
	// 		}
	// 	}

	// 	if(!goodNormal) {
	// 		throw MLException("Filter requires correct per vertex normals.<br>"
	// 							 "E.g. it is necessary that your <b>ALL</b> the input vertices have a proper, not-null normal.<br> "
	// 							 "Try enabling the <i>pre-clean<i> option and retry.<br><br>"
	// 							 "To permanently remove this problem:<br>"
	// 							 "If you encounter this error on a triangulated mesh try to use the <i>Remove Unreferenced Vertices</i> filter"
	// 							 "If you encounter this error on a pointcloud try to use the <i>Conditional Vertex Selection</i> filter"
	// 							 "with function '(nx==0.0) && (ny==0.0) && (nz==0.0)', and then <i>delete selected vertices</i>.<br>");
	// 	}

	// 	MeshModel *pm =md.addNewMesh("","Poisson mesh",false);
	// 	md.setVisible(pm->id(),false);
	// 	pm->updateDataMask(MeshModel::MM_VERTQUALITY);
	// 	if(goodColor)
	// 		pm->updateDataMask(MeshModel::MM_VERTCOLOR);

	// 	if(params.getBool("visibleLayer")) {
	// 		Box3m bb;
	// 		MeshModel *_mm=md.nextVisibleMesh();
	// 		while(_mm != nullptr){
	// 			bb.Add(_mm->cm.Tr,_mm->cm.bbox);
	// 			_mm=md.nextVisibleMesh(_mm);
	// 		}

	// 		MeshDocumentPointStream<Scalarm> documentStream(md);
	// 		_Execute<Scalarm,2,BOUNDARY_NEUMANN,PlyColorAndValueVertex<Scalarm> >(&documentStream,bb,pm->cm,pp,cb);
	// 	}
	// 	else {
	// 		MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
	// 		_Execute<Scalarm,2,BOUNDARY_NEUMANN,PlyColorAndValueVertex<Scalarm> >(&meshStream,md.mm()->cm.bbox,pm->cm,pp,cb);
	// 	}
	// 	pm->updateBoxAndNormals();
	// 	md.setVisible(pm->id(),true);
	// 	md.setCurrentMesh(pm->id());
	// 	if(currDirChanged)
	// 		QDir::setCurrent(currDir.path());
	// }
	// else {
	// 	wrongActionCalled(filter);
	// }
	return std::map<std::string, QVariant>();
}

RichParameterList FilterScreenedPoissonPlugin::initParameterList(
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
            parlist.addParam(RichInt ("pointWeight",
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

	}
	return parlist;
}

int FilterScreenedPoissonPlugin::postCondition(const QAction* filter) const
{
	if (ID(filter) == FP_SCREENED_POISSON){
		return MeshModel::MM_VERTNUMBER + MeshModel::MM_FACENUMBER;
	}
	else {
		return MeshModel::MM_ALL;
	}
}


FilterPlugin::FilterArity FilterScreenedPoissonPlugin::filterArity(const QAction*) const
{
	return VARIABLE;
}

