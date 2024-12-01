#define EIGEN_USE_GPU
#include <Eigen/Dense>
#include <cuda_runtime.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <vector>
#include<chrono>
using namespace Eigen;
using namespace std;

__global__
void Moving(int row,double* mp_gpu, double* mpd_gpu, double moving_point_step, double* ep_gpu, int* in_edge, 
    int* moving_point_cover_gpu, int* next_cover_gpu, int* before_cover_gpu, int* f, double* v,int rowv, 
    int* neighbor_gpu){
	Eigen::Map<Eigen::MatrixXd> moving_point(mp_gpu, row, 3);
	Eigen::Map<Eigen::MatrixXd> moving_point_direct(mpd_gpu, row, 3);
	Eigen::Map<Eigen::MatrixXd> edge_point(ep_gpu, row, 3);
    Eigen::Map<Eigen::MatrixXi> moving_point_cover(moving_point_cover_gpu, row, 1);
	Eigen::Map<Eigen::MatrixXi> before_cover(before_cover_gpu, row, 1);
	Eigen::Map<Eigen::MatrixXi> next_cover(next_cover_gpu, row, 1);
	Eigen::Map<Eigen::MatrixXi> F(f, row, 3);
	Eigen::Map<Eigen::MatrixXd> V(v, rowv, 3);
    Eigen::Map<Eigen::MatrixXi> neighbor(neighbor_gpu, row, 3);
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;
	for (int j = i;  j < row; j += stride) {
		if (j  < row) {
			if (moving_point.row(j) == RowVector3d(0, 0, 0)) {
				continue;
			}
			////判断是否到边
			if (in_edge[j] == 1) {
                int end = 0;
                in_edge[j] = 0;
                int cover = moving_point_cover(j, 0);
                int a = F(cover, 0);
                int b = F(cover, 1);
                int c = F(cover, 2);
                //BC
                double t1 = (moving_point(j, 1) - V(b, 1)) * (V(c, 2) - moving_point(j, 2)) - (moving_point(j, 2) - V(b, 2)) * (V(c, 1) - moving_point(j, 1));
                double t2 = (V(c, 1) - V(b, 1)) * moving_point_direct(j, 2) - (V(c, 2) - V(b, 2)) * moving_point_direct(j, 1);
                double t = t1 / t2;
                if (t > 0) {
                    RowVector3d point = moving_point.row(j) + t * moving_point_direct.row(j);
                    Eigen::RowVector3d v1 = V.row(b) - point;
                    Eigen::RowVector3d v2 = V.row(c) - point;

                    if (v1.dot(v2) < 0 && neighbor(cover, 2) != before_cover(j, 0)) {
                        before_cover(j, 0) = cover;
                        edge_point.row(j) = point;
                        next_cover(j, 0) = neighbor(cover, 2);
                        end = 1;
                    }

                }
                //AB
                t1 = (moving_point(j, 1) - V(a, 1)) * (V(b, 2) - moving_point(j, 2)) - (moving_point(j, 2) - V(a, 2)) * (V(b, 1) - moving_point(j, 1));
                t2 = (V(b, 1) - V(a, 1)) * moving_point_direct(j, 2) - (V(b, 2) - V(a, 2)) * moving_point_direct(j, 1);
                t = t1 / t2;
                if (t > 0 && end==0) {
                    RowVector3d point = moving_point.row(j) + t * moving_point_direct.row(j);
                    Eigen::RowVector3d v1 = V.row(a) - point;
                    Eigen::RowVector3d v2 = V.row(b) - point;
                    if (v1.dot(v2) < 0 && neighbor(cover, 0) != before_cover(j, 0)) {
                        before_cover(j, 0) = cover;
                        edge_point.row(j) = point;
                        next_cover(j, 0) = neighbor(cover, 0);
                        end = 1;
                    }

                }
                //AC
                t1 = (moving_point(j, 1) - V(a, 1)) * (V(c, 2) - moving_point(j, 2)) - (moving_point(j, 2) - V(a, 2)) * (V(c, 1) - moving_point(j, 1));
                t2 = (V(c, 1) - V(a, 1)) * moving_point_direct(j, 2) - (V(c, 2) - V(a, 2)) * moving_point_direct(j, 1);
                t = t1 / t2;
                if (t > 0 && end == 0) {
                    RowVector3d point = moving_point.row(j) + t * moving_point_direct.row(j);
                    Eigen::RowVector3d v1 = V.row(a) - point;
                    Eigen::RowVector3d v2 = V.row(c) - point;
                    if (v1.dot(v2) < 0 && neighbor(cover, 1) != before_cover(j, 0)) {
                        before_cover(j, 0) = cover;
                        edge_point.row(j) = point;
                        next_cover(j, 0) = neighbor(cover, 1);
                        end = 1;
                    }

                }
                if (end == 0) {
                    moving_point_direct.row(j) = RowVector3d(0, 0, 0);
                    moving_point.row(j) = RowVector3d(0, 0, 0);
                    edge_point.row(j) = RowVector3d(0, 0, 0);
                }
			}
			//计算现在的点的位置
			moving_point(j, 0) = moving_point(j, 0) + moving_point_direct(j, 0) * moving_point_step;
			moving_point(j, 1) = moving_point(j, 1) + moving_point_direct(j, 1) * moving_point_step;
			moving_point(j, 2) = moving_point(j, 2) + moving_point_direct(j, 2) * moving_point_step;
		}
	}
}


extern "C" void moving_point_gpu(int row, double* mp_gpu, double* mpd_gpu, double moving_point_step,
    double* ep_gpu, int* in_edge, int* moving_point_cover_gpu, int* next_cover_gpu, int* before_cover_gpu,
    int* f, double* v, int rowv, int* neighbor_gpu){
    Moving << <2, 256 >> > (row, mp_gpu, mpd_gpu, moving_point_step, 
        ep_gpu, in_edge, moving_point_cover_gpu, next_cover_gpu, 
        before_cover_gpu,f,v,rowv, neighbor_gpu);	
    cudaDeviceSynchronize();
	
}