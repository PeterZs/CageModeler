#pragma once

namespace Eigen {
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using MatrixXi = Matrix<int, -1, -1, 0, -1, -1>;
    using Vector3d = Matrix<double, 3, 1, 0, 3, 1>;
}

void calculatePCA(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& pca_basic_matrix, Eigen::MatrixXd& pca_based_mesh_vertices, Eigen::Vector3d& barycenter);

void computeBB(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const Eigen::MatrixXd& pca_basic_matrix, const Eigen::Vector3d& barycenter);

void voxelizeBB(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXd& cage_vertices, float degreeOfSparseness);

void generateCage(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces);

