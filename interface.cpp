#include"interface.h"

void vector_feild_visualization(Eigen::MatrixXd v, Eigen::MatrixXi f,Eigen::MatrixXd hh) {
    config Config;
    //string path = Config.path_of_mesh;
    v.transposeInPlace();
    f.transposeInPlace();
    graph G(v,f,hh);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.core().camera_zoom = Config.camera_zoom;
    G.show_graph(viewer);
    G.show_grad(viewer);
    cout << "11111111111" << endl;
    G.show_point(viewer);

    // Define timing movements and resets
    double last_time = glfwGetTime();
    double start_tiem = glfwGetTime();
    //Modify render operation
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& viewer) {
        viewer.core().is_animating = true;
        double t = glfwGetTime();
        if (t - start_tiem >= Config.time_restart) {
            viewer.data().clear_points();
            G.restart();
            start_tiem = glfwGetTime();
            return false;
        }
        if (t - last_time >= Config.time_per_step) {
            G.check_point_in_edge();
            G.move_point(viewer);
            last_time = t;
        }
        return false;
        };

    viewer.launch();
}
