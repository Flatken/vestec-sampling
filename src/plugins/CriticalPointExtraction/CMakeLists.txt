vestec_add_paraview_plugin(vestecCriticalPointExtraction
	SOURCES VestecCriticalPointExtractionAlgorithm.cxx
	PLUGIN_VERSION 1.0
	PLUGIN_XML VestecCriticalPointExtraction.xml
	LINK Eigen3::Eigen Cnl::Cnl)
