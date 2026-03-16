#pragma once

#include <Eigen/Geometry>
#include <vector>

struct BBVoxel {
    Eigen::Vector3d start_pt;
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> centers, voxel_pts;
    int n_voxel[3];
    float res_voxel[3];
    std::vector<std::vector<std::vector<int8_t>>> voxel_types; // 0 = feature voxel, -1 = inner voxel, 1 = outer voxel
    std::vector<int> non_manifold_vertices;
    std::vector<Eigen::Vector3d> splitted_vertices;
};

// struct Hit
// {
//     double t;                  
//     Eigen::Vector3d position;  
//     Eigen::Vector3d normal;    
//     int face_id;               
// };


void calculatePCA(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& pca_basic_matrix, Eigen::MatrixXd& pca_based_mesh_vertices, Eigen::Vector3d& barycenter);

void computeBB(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const Eigen::MatrixXd& pca_basic_matrix, const Eigen::Vector3d& barycenter);

void voxelizeBB(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXd& cage_vertices, BBVoxel& voxels, float degreeOfSparseness);

void minMax(const double v0, const double v1, const double v2, double& min, double& max);
bool planeBoxOverlap(const Eigen::Vector3d normal, const double d, const Eigen::Vector3d maxbox);
bool axistest_x01(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
bool axistest_y02(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
bool axistest_z12(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
bool axistest_z0(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
bool axistest_x2(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
bool axistest_y1(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox);
void identifyFeatureVoxels(BBVoxel& voxels, const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXi& mesh_faces);

void identifyOuterVoxels(BBVoxel& voxels);
void identifyOuterVoxelsWithFillingAlgo(BBVoxel& voxels, const Eigen::Vector3i seed);

void extractOuterSurface(BBVoxel& voxels, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces);

void smoothCage(Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, float lambda_smooth, int it_smooth, Eigen::MatrixXd& pca_based_mesh_vertices, const Eigen::MatrixXi& mesh_faces);

void generateCageCoarseBouding(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXi& mesh_faces, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const float sparse_factor);

