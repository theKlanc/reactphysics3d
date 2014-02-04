# Minimum cmake version required
cmake_minimum_required(VERSION 2.6)

# Project configuration
PROJECT(Fluid)

# Where to build the executables
SET(EXECUTABLE_OUTPUT_PATH "${OUR_EXECUTABLE_OUTPUT_PATH}/fluid")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${EXECUTABLE_OUTPUT_PATH})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${EXECUTABLE_OUTPUT_PATH})

# Copy the shaders used for the demo into the build directory
FILE(COPY "${OPENGLFRAMEWORK_DIR}/src/shaders/" DESTINATION "${EXECUTABLE_OUTPUT_PATH}/shaders/")

# Copy the meshes used for the demo into the build directory
FILE(COPY "../common/meshes/" DESTINATION "${EXECUTABLE_OUTPUT_PATH}/meshes/")

# Headers
INCLUDE_DIRECTORIES("${OPENGLFRAMEWORK_DIR}/src/" "../common/glfw/include/" "../common/")

# Source files
SET(FLUID_SOURCES
    Fluid.cpp
    Scene.cpp
    Scene.cpp
    Scene.h
    "../common/Sphere.cpp"
    "../common/Sphere.h"
    "../common/Box.cpp"
    "../common/Box.h"
    "../common/Viewer.cpp"
    "../common/Viewer.h"
)

# Create the executable
ADD_EXECUTABLE(fluid ${FLUID_SOURCES})

# Link with libraries
TARGET_LINK_LIBRARIES(fluid reactphysics3d openglframework glfw ${GLFW_LIBRARIES})
