#include <cagegenerations/CoarseBoundingCagesByXian.h>

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <limits>

static const double EPSILON_BB = 0.15;


void calculatePCA(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& pca_basic_matrix, Eigen::MatrixXd& pca_based_mesh_vertices, Eigen::Vector3d& barycenter)
{
    barycenter = mesh_vertices.colwise().mean(); // barycenter of mesh

    Eigen::MatrixXd X = mesh_vertices.rowwise() - barycenter.transpose();

    int num_vertices = mesh_vertices.rows();

    Eigen::Matrix3d covariance_mat = (X.transpose() * X) / (num_vertices - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance_mat); // eigendecomposition
    Eigen::Matrix3d eigenvectors = es.eigenvectors();
    // Eigen::Vector3d eigenvalues = es.eigenvalues();

    pca_basic_matrix = Eigen::MatrixXd::Zero(3,3);
    for (int i = 0; i < eigenvectors.cols(); i++) 
    {   
        Eigen::Vector3d vec = eigenvectors.col(i).real();
        pca_basic_matrix.col(i) = vec / vec.norm();
    }
    
    pca_based_mesh_vertices = mesh_vertices * pca_basic_matrix;
}

void computeBB(const Eigen::MatrixXd& vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const Eigen::MatrixXd& pca_basic_matrix, const Eigen::Vector3d& barycenter)
{
    // calculate max and min coordinates for BB
    Eigen::Vector3d max_coordinates(-std::numeric_limits<double>::max(),   -std::numeric_limits<double>::max(),    -std::numeric_limits<double>::max());
    Eigen::Vector3d min_coordinates(std::numeric_limits<double>::max(),    std::numeric_limits<double>::max(),     std::numeric_limits<double>::max());
    Eigen::Vector3d p;
    for (int i = 0; i < vertices.rows(); i++) 
    {   
        p = vertices.row(i);
        for (int j = 0; j < 3; j++) 
        {
            if (p[j] > max_coordinates[j])  max_coordinates[j] = p[j];
            if (p[j] < min_coordinates[j])  min_coordinates[j] = p[j];
        }
    }
    for (int j = 0; j < 3; j++) 
    {   
        double d = max_coordinates[j] - min_coordinates[j];
        max_coordinates[j] += EPSILON_BB * d;
        min_coordinates[j] -= EPSILON_BB * d;
    }

    // calculate actual BB
    // and store it already in the cage for instance
    // delete cage for now
    cage_vertices = Eigen::MatrixXd::Zero(8, 3);
    cage_faces = Eigen::MatrixXi::Zero(12, 8);

    cage_vertices.row(0) = Eigen::Vector3d(min_coordinates[0], min_coordinates[1], min_coordinates[2]);
    cage_vertices.row(1) = Eigen::Vector3d(max_coordinates[0], min_coordinates[1], min_coordinates[2]);
    cage_vertices.row(2) = Eigen::Vector3d(max_coordinates[0], min_coordinates[1], max_coordinates[2]);
    cage_vertices.row(3) = Eigen::Vector3d(min_coordinates[0], min_coordinates[1], max_coordinates[2]);
    cage_vertices.row(4) = Eigen::Vector3d(min_coordinates[0], max_coordinates[1], min_coordinates[2]);
    cage_vertices.row(5) = Eigen::Vector3d(max_coordinates[0], max_coordinates[1], min_coordinates[2]);
    cage_vertices.row(6) = Eigen::Vector3d(max_coordinates[0], max_coordinates[1], max_coordinates[2]);
    cage_vertices.row(7) = Eigen::Vector3d(min_coordinates[0], max_coordinates[1], max_coordinates[2]);

    cage_faces.row(0) = Eigen::Vector3i(0,1,2);
    cage_faces.row(1) = Eigen::Vector3i(0,2,3);
    cage_faces.row(2) = Eigen::Vector3i(1,5,6);
    cage_faces.row(3) = Eigen::Vector3i(1,6,2);
    cage_faces.row(4) = Eigen::Vector3i(5,4,7);
    cage_faces.row(5) = Eigen::Vector3i(5,7,6);
    cage_faces.row(6) = Eigen::Vector3i(4,0,3);
    cage_faces.row(7) = Eigen::Vector3i(4,3,7);
    cage_faces.row(8) = Eigen::Vector3i(4,5,1);
    cage_faces.row(9) = Eigen::Vector3i(4,1,0);
    cage_faces.row(10) = Eigen::Vector3i(3,2,6);
    cage_faces.row(11) = Eigen::Vector3i(3,6,7);

    // transform back
    cage_vertices = cage_vertices * pca_basic_matrix.transpose();
    for (int i = 0; i < cage_vertices.rows(); i++)
    {   
        cage_vertices.row(i) += barycenter;
    }
}

void generateCage(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{   
    
    Eigen::MatrixXd pca_basic_matrix;
    Eigen::MatrixXd pca_based_mesh_vertices;
    Eigen::Vector3d barycenter;
    calculatePCA(mesh_vertices, pca_basic_matrix, pca_based_mesh_vertices, barycenter);

    computeBB(pca_based_mesh_vertices, cage_vertices, cage_faces, pca_basic_matrix, barycenter);

    
};