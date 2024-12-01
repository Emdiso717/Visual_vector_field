cmake_minimum_required(VERSION 3.16)
project(project)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
# Other modules you could enable
#igl_include(embree)
#igl_include(opengl)
#igl_include(stb)
#igl_include(predicates)
#igl_include(xml)
#igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
#igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

enable_language(CUDA)

# Add your project files
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/data)
include_directories(${PROJECT_SOURCE_DIR}/eigen-3.4.0)
file(GLOB SRC_FILES *.cpp)
file(GLOB CUDA_SRC_FILES *.cu)
add_executable(${PROJECT_NAME} "${SRC_FILES}"
"${CUDA_SRC_FILES}")
# Link igl (and the glfw module) to your project
target_link_libraries(${PROJECT_NAME} PUBLIC 
  igl::glfw    igl::imgui
  ## Other modules you could link to
  # igl::embree
  # igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  # igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  # igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  # igl_restricted::triangle
  )

# Find CUDA libraries
find_package(CUDA REQUIRED)

# Include CUDA Toolkit include directories
include_directories(${CUDAToolkit_INCLUDE_DIRS})

# Link CUDA libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${CUDA_LIBRARIES})
