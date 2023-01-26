# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES	filter_iPSR.cpp)

set(HEADERS
	src/Allocator.h
	src/Array.h
	src/BSplineData.h
	src/BinaryNode.h
	src/CmdLineParser.h
	src/FEMTree.h
	src/Factor.h
	src/FunctionData.h
	src/Geometry.h
	src/LinearSolvers.h
	src/MAT.h
	src/MarchingCubes.h
	src/MyMiscellany.h
	src/PPolynomial.h
	src/Ply.h
	src/PointStream.h
	src/PointStreamData.h
	src/PoissonRecon.h
	src/Polynomial.h
	src/RegularTree.h
	src/SparseMatrix.h
	src/SparseMatrixInterface.h
	src/Window.h
	src/kdtree.h
	src/utility.h
	filter_iPSR.h
	#poisson_utils.h
	)

set(sNL_HEADERS
	src/Array.inl
	src/BMPStream.inl
	src/BSplineData.inl
	src/CmdLineParser.inl
	src/FEMTree.Evaluation.inl
	src/FEMTree.Initialize.inl
	src/FEMTree.IsoSurface.specialized.inl
	src/FEMTree.SortedTreeNodes.inl
	src/FEMTree.System.inl
	src/FEMTree.WeightedSamples.inl
	src/FEMTree.inl
	src/FunctionData.inl
	src/Geometry.inl
	src/MAT.inl
	src/PPolynomial.inl
	src/PointStream.inl
	src/PoissonRecon.inl
	src/Polynomial.inl
	src/RegularTree.inl
	src/SparseMatrix.inl
	src/SparseMatrixInterface.inl
	src/Window.inl)
add_meshlab_plugin(filter_iPSR ${SOURCES} ${HEADERS} ${INL_HEADERS})

target_compile_definitions(filter_iPSR PRIVATE BRUNO_LEVY_FIX
														   FOR_RELEASE)

if(OpenMP_CXX_FOUND)
	target_link_libraries(filter_iPSR PRIVATE OpenMP::OpenMP_CXX)
endif()
