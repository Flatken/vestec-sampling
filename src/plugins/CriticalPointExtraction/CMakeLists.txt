vestec_add_paraview_plugin(vestecCriticalPointExtraction
	SOURCES VestecCriticalPointExtractionAlgorithm.cxx
	HEADERS VestecCriticalPointExtractionAlgorithm.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML VestecCriticalPointExtraction.xml
	LINKS Eigen3::Eigen OpenMP::OpenMP_CXX vtkm_cont)
