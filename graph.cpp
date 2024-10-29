
#include "graph.h"

void graph::get_neighbor()
{
    neighbor = Eigen::MatrixXi(F.rows(), 3);
    string file_name = B.path_of_mesh + ".txt";
    std::ifstream file(file_name);
    string line = "\0";
    getline(file, line);
    if (line == B.path_of_mesh) {
        int i = 0;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double value;
            int j = 0;
            while (iss >> value) {
                neighbor(i, j) = value;
                j++;
            }
            i++;
        }
        file.close();
    }else{
        file.close();
        ofstream file1(file_name);
        file1 << B.path_of_mesh << endl;;
        // [0,1] [0,2] [1,2]
        for (int i = 0; i < F.rows(); i++) {
            for (int j = 0; j < F.rows(); j++) {
                //除了本面
                if (j == i)
                    continue;
                if (F(j, 0) == F(i, 0)) {
                    if (F(j, 1) == F(i, 1) || F(j, 2) == F(i, 1))
                        neighbor(i, 0) = j;
                    else if (F(j, 1) == F(i, 2) || F(j, 2) == F(i, 2))
                        neighbor(i, 1) = j;
                }
                else if (F(j, 1) == F(i, 0)) {
                    if (F(j, 0) == F(i, 1) || F(j, 2) == F(i, 1))
                        neighbor(i, 0) = j;
                    else if (F(j, 0) == F(i, 2) || F(j, 2) == F(i, 2))
                        neighbor(i, 1) = j;
                }
                else if (F(j, 2) == F(i, 0)) {
                    if (F(j, 0) == F(i, 1) || F(j, 1) == F(i, 1))
                        neighbor(i, 0) = j;
                    else if (F(j, 0) == F(i, 2) || F(j, 1) == F(i, 2))
                        neighbor(i, 1) = j;
                }
                else {
                    if (F(j, 0) == F(i, 1)) {
                        if (F(j, 1) == F(i, 2) || F(j, 2) == F(i, 2))
                            neighbor(i, 2) = j;
                    }
                    else if (F(j, 1) == F(i, 1)) {
                        if (F(j, 0) == F(i, 2) || F(j, 2) == F(i, 2))
                            neighbor(i, 2) = j;
                    }
                    else if (F(j, 2) == F(i, 1)) {
                        if (F(j, 0) == F(i, 2) || F(j, 1) == F(i, 2))
                            neighbor(i, 2) = j;
                    }
                }
            }
            file1 << neighbor(i, 0) << ' ' << neighbor(i, 1) << ' ' << neighbor(i, 2) << endl;
        }
        file.close();
    }
}


double graph::cal_next_edge_point(int i, int a, int b)
{
    double t1 = (moving_point(i, 1) - V(a, 1)) * (V(b, 2) - moving_point(i, 2)) - (moving_point(i, 2) - V(a, 2)) * (V(b, 1) - moving_point(i, 1));
    double t2 = (V(b, 1) - V(a, 1)) * moving_point_direct(i, 2) - (V(b, 2) - V(a, 2)) * moving_point_direct(i, 1);
    double t3 = t1 / t2;
    return t3;
}
void graph::add_point_in_empty_cover()
{
    int a1 = cover_without_point.size();
    int a2 = point_deleted.size();
    int num_add = min(a1, a2);
    for (int i = 0; i < num_add; i++) {
        int new_point_number = point_deleted[i];
        int new_cover_number = cover_without_point[i];
        point_in_cover[new_cover_number].push_back(new_point_number);
        moving_point_direct.row(new_point_number) = K.row(new_cover_number);
        moving_point.row(new_point_number) = C.row(new_cover_number);
        moving_point_cover(new_point_number, 0) = new_cover_number;
        before_cover(new_point_number, 0) = -1;
        cal_edge_point(new_point_number);
    }
    point_deleted.erase(point_deleted.begin(), point_deleted.begin()+ num_add);
    cover_without_point.erase(cover_without_point.begin(), cover_without_point.begin()+ num_add);
}
graph::graph(string path)
{   
    //read from file
    if (path.find(".obj") != std::string::npos) 
    { 
        igl::readOBJ(path, V, F); 
    }
    else {
        igl::readOFF(path, V, F);
    }
    H = Eigen::MatrixXd(V.rows(), 1);
    std::ifstream file(B.path_of_grad);
    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
        H(i,0) = std::stod(line);
        i++;
    }
    // cal grad
    K = Eigen::MatrixXd(F.rows(), 3);
    grad();
    // find centers
    C = Eigen::MatrixXd(F.rows(), 3);
    next_cover = Eigen::MatrixXi(F.rows(), 1);
    igl::barycenter(V, F, C);
    //basic of moving points
    moving_point_direct = K;
    moving_point = C;
    moving_point_cover = Eigen::MatrixXi(F.rows(), 1);
    before_cover = Eigen::MatrixXi(F.rows(), 1);

    //init of points in a cover
    point_in_cover = new vector<int>[F.rows()];


    //find neighbors
    graph::get_neighbor();
    //find the first edge points
    edge_point = MatrixXd(F.rows(), 3);
    for (int i = 0; i < F.rows(); i++) {
        moving_point_cover(i, 0) = i;
        before_cover(i, 0) = -1;
        //init the first cover condition
        (point_in_cover[i]).push_back(i);
    }
    for (int i = 0; i < F.rows(); i++) {
        graph::cal_edge_point(i);
    }
}



void graph::grad()
{
    for (int i = 0; i < F.rows(); i++)
    {
        // init some varible
        int a1 = F(i, 0);
        int a2 = F(i, 1);
        int a3 = F(i, 2);
        RowVector3d v1 = V.row(a1) - V.row(a2);
        RowVector3d v2 = V.row(a3) - V.row(a1);
        RowVector3d v3 = V.row(a2) - V.row(a3);
        // calculate normal vector
        RowVector3d n = v1.cross(v2);
        double s = 0.5 * n.norm();
        n = n / n.norm();
        // a1
        RowVector3d va1 = -v3.cross(n);
        va1 = (va1 / va1.norm()) / (s / v3.norm());
        // a2
        RowVector3d va2 = -v2.cross(n);
        va2 = (va2 / va2.norm()) / (s / v2.norm());
        // a3
        RowVector3d va3 = -v1.cross(n);
        va3 = (va3 / va3.norm()) / (s / v1.norm());
        RowVector3d grad = va1 * H(a1, 0) + va2 * H(a2, 0) + va3 * H(a3, 0);
        grad = grad / grad.norm();
        K.row(i) = grad;
    }
}
void graph::show_graph(igl::opengl::glfw::Viewer& viewer)
{
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
}
void graph::show_grad(igl::opengl::glfw::Viewer& viewer)
{
    viewer.data().set_data(H);
   // int num_colors = 256;
   // Eigen::MatrixXd custom_colormap(num_colors, 3);
   // for (int i = 0; i < num_colors; ++i) {
   //     float t = static_cast<float>(i) / (num_colors - 1); 
   //     custom_colormap(i, 0) = 1.0f - t; 
   //     custom_colormap(i, 1) = 1.0f - t; 
   //     custom_colormap(i, 2) = 2.0f * t; 
   // }
   // viewer.data().set_colormap(custom_colormap);
   viewer.data().add_edges(C, C + B.grad_line_length * K, B.grad_line_color);

}
void graph::show_point(igl::opengl::glfw::Viewer& viewer)
{
    viewer.data().point_size = B.point_size;
    double dis_each_point = B.length_of_point/B.num_points;
    float p_color = 1.0 / B.num_points;
    for (int i = 0; i < B.num_points; i++) {
        viewer.data().add_points(moving_point - dis_each_point * moving_point_direct * i, B.point_color + Eigen::RowVector3d(0,p_color,p_color)*i);
    }
}
void graph::move_point(igl::opengl::glfw::Viewer& viewer)
{
    moving_point = moving_point + moving_point_direct * B.moving_point_step;
    viewer.data().clear_points();
    graph::show_point(viewer);
}


void graph::cal_edge_point(int i)
{
    if (moving_point_direct(i,0)==0 && moving_point_direct(i, 1) == 0 && moving_point_direct(i, 2) == 0)
        return;
    //确定平面
    int cover = moving_point_cover(i, 0);
    int a = F(cover, 0);
    int b = F(cover, 1);
    int c = F(cover, 2);
    //BC
    double t = cal_next_edge_point(i, b, c);
    if (t > 0) {
        RowVector3d point = moving_point.row(i) + t * moving_point_direct.row(i);
        Eigen::RowVector3d v1 = V.row(b) - point;
        Eigen::RowVector3d v2 = V.row(c) - point;
        
        if (v1.dot(v2) < 0 && neighbor(cover, 2) != before_cover(i,0)) {
            before_cover(i, 0) = cover;
            edge_point.row(i) = point;
            next_cover(i, 0) = neighbor(cover,2);
            return;
        }

    }
    //AB
    t = cal_next_edge_point(i, a, b);
    if (t > 0) {
        RowVector3d point = moving_point.row(i) + t * moving_point_direct.row(i);
        Eigen::RowVector3d v1 = V.row(a) - point;
        Eigen::RowVector3d v2 = V.row(b) - point;
        if (v1.dot(v2) < 0 && neighbor(cover, 0) != before_cover(i, 0)) {
            before_cover(i, 0) = cover;
            edge_point.row(i) = point;
            next_cover(i, 0) = neighbor(cover, 0);
            return;
        }

    }
    //AC
    t = cal_next_edge_point(i, a, c);
    if (t > 0) {
        RowVector3d point = moving_point.row(i) + t * moving_point_direct.row(i);
        Eigen::RowVector3d v1 = V.row(a) - point;
        Eigen::RowVector3d v2 = V.row(c) - point;
        if (v1.dot(v2) < 0 && neighbor(cover, 1) != before_cover(i, 0)) {
            before_cover(i, 0) = cover;
            edge_point.row(i) = point;
            next_cover(i, 0) = neighbor(cover, 1);
            return;
        }

    }
    auto newEnd = std::remove(point_in_cover[moving_point_cover(i, 0)].begin(), point_in_cover[moving_point_cover(i, 0)].end(), i);
    point_in_cover[moving_point_cover(i, 0)].erase(newEnd, point_in_cover[moving_point_cover(i, 0)].end());
    if (point_in_cover[moving_point_cover(i, 0)].size() == 0) {
        auto temp_element = find(cover_without_point.begin(), cover_without_point.end(), moving_point_cover(i, 0));
        if (temp_element == cover_without_point.end()) {
            cover_without_point.push_back(moving_point_cover(i, 0));
        }
    }
    //高矢量边
    moving_point_direct.row(i) = RowVector3d(0, 0, 0);
    moving_point.row(i) = RowVector3d(0, 0, 0);
    edge_point.row(i) = RowVector3d(0, 0, 0);
    point_deleted.push_back(i);
}

void graph::check_point_in_edge()
{
    for (int i = 0; i < F.rows(); i++) {
        if (moving_point.row(i) == RowVector3d(0, 0, 0)) {
            continue;
        }
        RowVector3d v1 = moving_point.row(i);
        RowVector3d v2 = edge_point.row(i);
        if ((v1 - v2).norm() <= B.moving_point_step) {
            //delete from the former cover
            auto newEnd = std::remove(point_in_cover[moving_point_cover(i, 0)].begin(), point_in_cover[moving_point_cover(i, 0)].end(), i);
            (point_in_cover[moving_point_cover(i, 0)]).erase(newEnd, point_in_cover[moving_point_cover(i, 0)].end());
            //add to the new cover
            (point_in_cover[next_cover(i, 0)]).push_back(i);
            auto find_element = find(cover_without_point.begin(), cover_without_point.end(), next_cover(i, 0));
            if (find_element != cover_without_point.end()) {
                cover_without_point.erase(find_element);
            }
            if (point_in_cover[moving_point_cover(i, 0)].size() == 0) {
                auto temp_element = find(cover_without_point.begin(), cover_without_point.end(), moving_point_cover(i, 0));
                if (temp_element == cover_without_point.end()) {
                    cover_without_point.push_back(moving_point_cover(i, 0));
                }
            }
            if (point_in_cover[next_cover(i, 0)].size() >= 5) {
                auto newEnd = std::remove(point_in_cover[next_cover(i, 0)].begin(), point_in_cover[next_cover(i, 0)].end(), i);
                (point_in_cover[next_cover(i, 0)]).erase(newEnd, point_in_cover[next_cover(i, 0)].end());
                moving_point_direct.row(i) = RowVector3d(0, 0, 0);
                moving_point.row(i) = RowVector3d(0, 0, 0);
                edge_point.row(i) = RowVector3d(0, 0, 0);
                point_deleted.push_back(i);
                continue;
            }
            moving_point.row(i) = edge_point.row(i);
            moving_point_cover(i, 0) = next_cover(i, 0);
            moving_point_direct.row(i) = K.row(moving_point_cover(i, 0));
            cal_edge_point(i);
        }
    }
    add_point_in_empty_cover();
}

void graph::restart()
{
    moving_point = C;
    moving_point_direct = K;
    point_deleted.clear();
    cover_without_point.clear();
    for (int i = 0; i < F.rows(); i++) {
        moving_point_cover(i, 0) = i;
        before_cover(i, 0) = -1;
        point_in_cover[i].clear();
        point_in_cover[i].push_back(i);
    }
    for (int i = 0; i < F.rows(); i++) {
        graph::cal_edge_point(i);
    }
}

void graph::debug(igl::opengl::glfw::Viewer& viewer)
{
    viewer.data().point_size = B.point_size;
    viewer.data().add_points(C, RowVector3d(0, 1, 0));
    //viewer.data().add_points(edge_point, RowVector3d(1, 0, 0));
    moving_point = edge_point;
    moving_point_cover = next_cover;
    for (int i = 0; i < F.rows(); i++) {
        moving_point_direct.row(i) = K.row(next_cover(i, 0));
        viewer.data().add_edges(edge_point.row(i), edge_point.row(i) + B.grad_line_length * K.row(next_cover(i, 0)), RowVector3d(0, 0, 0));
    }
    for (int i = 0; i < F.rows(); i++) {
        cal_edge_point(i);
    }
    viewer.data().add_points(edge_point, RowVector3d(0, 0, 1));

}


