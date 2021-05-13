cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
set(OWN_GXX_FLAGS "-Wextra -fstack-protector-all -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS} -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
	utils/easy_image.h
	utils/easy_image.cc
	utils/ini_configuration.h
	utils/ini_configuration.cc
	utils/utils.h
	utils/utils.cpp
	utils/l_parser.h
	utils/l_parser.cc

	vector/vector3d.h
	vector/vector3d.cc

		Geometry/L2D.h
		Geometry/L2D.cpp
		Geometry/L3D.h
		Geometry/L3D.cpp
        ImageGeneration/ImageSpecifications.h
        ImageGeneration/ImageSpecifications.cpp
        ImageGeneration/ImageGenerator.cpp
        ImageGeneration/ImageGenerator.h

	engine.cc)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/build )
