#include "main.h"



int main(int argc, char *argv[])
{
    // Read from OBJ file just for test
    ifstream file("../input.json");
    string json_str((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    string root = "../data/";
    file.close();
    json json_obj = json::parse(json_str);
    string path_of_mesh = json_obj["path_of_mesh"];
    string path_of_grad = json_obj["path_of_grad"];
    string path = root + path_of_mesh;
    path_of_grad = root + path_of_grad;
    Eigen::MatrixXd v;
    Eigen::MatrixXi f;
    if (path.find(".obj") != std::string::npos) 
    { 
        igl::readOBJ(path, v, f); 
    }
    else {
        igl::readOFF(path, v, f);
    }
    Eigen::MatrixXd hh = Eigen::MatrixXd(v.rows(), 1);
    std::ifstream file1(path_of_grad);
    std::string line;
    int i = 0;
    while (std::getline(file1, line)) {
        hh(i, 0) = std::stod(line);
        i++;
    }
    f.transposeInPlace();
    v.transposeInPlace();
	vector_feild_visualization(v,f,hh);
}
