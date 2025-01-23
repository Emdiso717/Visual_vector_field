#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/barycenter.h>
#include <iostream>
#include <igl/readOFF.h>
#include <imgui.h>
#include <Eigen/Core>
#include<vector>
#include <igl/ray_mesh_intersect.h>
#include "graph.h"
#include <nlohmann/json.hpp>
#include <cuda_runtime.h>
#include "interface.h"


using namespace Eigen;
using namespace std;

