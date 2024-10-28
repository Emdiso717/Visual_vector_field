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

using namespace Eigen;
using namespace std;

struct config
{
    string path_of_mesh;
    double time_per_step;
    double camera_zoom;
    double time_restart;
    config() {
        ifstream file("../input.json");
        string json_str((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
        string root = "../data/";
        file.close();
        json json_obj = json::parse(json_str);
        path_of_mesh = json_obj["path_of_mesh"];
        path_of_mesh = root + path_of_mesh;
        time_per_step = json_obj["time_per_step"];
        camera_zoom = json_obj["camera_zoom"];
        time_restart = json_obj["time_restart"];
    }
};