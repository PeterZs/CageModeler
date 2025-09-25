#pragma once

#include <Eigen/Geometry>
#include <vector>

struct BBVoxel {
    Eigen::Vector3d start_pt;
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> centers, voxel_pts;
    int n_voxel[3];
    float res_voxel[3];
};

void calculatePCA(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& pca_basic_matrix, Eigen::MatrixXd& pca_based_mesh_vertices, Eigen::Vector3d& barycenter);

void computeBB(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const Eigen::MatrixXd& pca_basic_matrix, const Eigen::Vector3d& barycenter);

void voxelizeBB(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXd& cage_vertices, BBVoxel& voxels, float degreeOfSparseness);

void generateCage(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces);

