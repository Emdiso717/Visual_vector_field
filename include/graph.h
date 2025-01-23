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
#include <vector>
#include <igl/ray_mesh_intersect.h>
#include <nlohmann/json.hpp>
#include <random>
#include <cuda_runtime.h>
#include <chrono>


using json = nlohmann::json;
using namespace Eigen;
using namespace std;

struct basic
{
	RowVector3d grad_line_color;
	float grad_line_length;
	RowVector3d point_color;
	int num_points;
	double point_size;
	double length_of_point;
	double moving_point_step;
	string path_of_mesh;
	string path_of_grad;
	int show_grad_line;
	basic() {
		ifstream file("../input.json");
		string root = "../data/";
		string json_str((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
		file.close();
		json json_obj = json::parse(json_str);
		grad_line_length = json_obj["grad_line_length"];
		num_points = json_obj["num_points"];
		point_size = json_obj["point_size"];
		length_of_point = json_obj["length_of_point"];
		grad_line_color = RowVector3d(json_obj["grad_line_color"][0], json_obj["grad_line_color"][1], json_obj["grad_line_color"][2]);
		point_color = RowVector3d(json_obj["point_color"][0], json_obj["point_color"][1], json_obj["point_color"][2]);
		moving_point_step= json_obj["moving_point_step"];
		path_of_mesh = json_obj["path_of_mesh"];
		path_of_grad = json_obj["path_of_grad"];
		show_grad_line = json_obj["show_grad_line"];
		path_of_mesh = root + path_of_mesh;
		path_of_grad = root + path_of_grad;
		cout << "Files that need to be rendered:" << path_of_mesh << endl;
	}
};


class graph
{
	private:
		/*ͼ������ �������ı�*/
		Eigen::MatrixXd V;	//��
		Eigen::MatrixXi F;	//��
		Eigen::MatrixXd K;	//������
		Eigen::MatrixXd C;	//�е�
		Eigen::MatrixXd H;	//��ʸ��ֵ
		Eigen::MatrixXi neighbor;	//���ھ�
		Eigen::MatrixXd change_cover;
		int* neighbor_gpu;
		int* f;
		double * v;

		/*�˶�������� �仯*/
		int* next_cover_gpu;	//���˶���һ����
		int* before_cover_gpu, * moving_point_cover_gpu;	//���˶�ǰһ���棬���ڵ���
		double* mp_gpu, * mpd_gpu;	//������꣬���˶�����
		double* ep_gpu ;	//������һ���ཻ�ߵĽ���
		int *in_edge;	//�ж��Ƿ񵽱�
	

		/*��¼ÿ�����еĵ�*/
		vector<int> *point_in_cover;
		vector<int> point_deleted;
		vector<int> cover_without_point;

		basic B;
		void get_neighbor();
		double cal_next_edge_point(int i, int a, int b);
		void add_point_in_empty_cover();
	public:
		graph(Eigen::MatrixXd vv, Eigen::MatrixXi ff,Eigen::MatrixXd hh);
		void grad();
		void show_graph(igl::opengl::glfw::Viewer& viewer);
		void show_grad(igl::opengl::glfw::Viewer& viewer);
		void show_point(igl::opengl::glfw::Viewer& viewer);
		void move_point(igl::opengl::glfw::Viewer& viewer);
		void cal_edge_point(int i);
		void check_point_in_edge();
		void restart();


};

