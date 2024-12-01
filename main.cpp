#include "main.h"



int main(int argc, char *argv[])
{
   config Config;
   string path = Config.path_of_mesh;
   graph G(path);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.core().camera_zoom = Config.camera_zoom;
  G.show_graph(viewer);
  G.show_grad(viewer);
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
          G.move_point(viewer);
          G.check_point_in_edge();
          last_time = t;
      }
      return false; 
   };
  viewer.launch();
}
