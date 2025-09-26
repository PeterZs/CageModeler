#include <cagegenerations/CoarseBoundingCagesByXian.h>

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <limits>

static const double EPSILON_BB = 0.15;
static const float PREDEFINED_SPARSE_FACTOR = 0.5;


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
    cage_faces = Eigen::MatrixXi::Zero(12, 3);

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

void voxelizeBB(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXd& cage_vertices, BBVoxel& voxels, float degreeSparseness = PREDEFINED_SPARSE_FACTOR)
{    
    Eigen::Vector3d max_coordinates(-std::numeric_limits<double>::max(),   -std::numeric_limits<double>::max(),    -std::numeric_limits<double>::max());
    Eigen::Vector3d min_coordinates(std::numeric_limits<double>::max(),    std::numeric_limits<double>::max(),     std::numeric_limits<double>::max());
    Eigen::Vector3d p;
    for (int i = 0; i < cage_vertices.rows(); i++) 
    {   
        p = cage_vertices.row(i);
        for (int j = 0; j < 3; j++) 
        {
            if (p[j] > max_coordinates[j])  max_coordinates[j] = p[j];
            if (p[j] < min_coordinates[j])  min_coordinates[j] = p[j];
        }
    }
    int largest_dim = 0;
    float largest_length = max_coordinates[0] - min_coordinates[0];
    for (int i = 1; i < 3; i++) 
    {
        float diff = max_coordinates[i] - min_coordinates[i];
        if (largest_length < diff)
        {
            largest_dim = i;
            largest_length = diff;
        }
    }

    int num_of_voxel_each_dim[3];
    float res_of_voxel_each_dim[3];
    num_of_voxel_each_dim[largest_dim] = (int)(sqrt(mesh_vertices.rows() * degreeSparseness) / 6.0);
    res_of_voxel_each_dim[largest_dim] = (max_coordinates[largest_dim] - min_coordinates[largest_dim]) / num_of_voxel_each_dim[largest_dim];
    for (int i = 0; i < 3; i++)
    {
        if (i != largest_dim) 
        {
            float diff = max_coordinates[i] - min_coordinates[i];
            num_of_voxel_each_dim[i] = (int) diff / res_of_voxel_each_dim[largest_dim] + 2;
            res_of_voxel_each_dim[i] = diff / num_of_voxel_each_dim[i];
        }
    }

    // Calculation of the center of all voxels and the voxel grid points
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> centers;
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> voxel_pts;
    std::vector<std::vector<std::vector<int8_t>>> voxel_types;
    voxels = BBVoxel{
        min_coordinates,
        centers,
        voxel_pts,
        {num_of_voxel_each_dim[0], num_of_voxel_each_dim[1], num_of_voxel_each_dim[2]},
        {res_of_voxel_each_dim[0], res_of_voxel_each_dim[1], res_of_voxel_each_dim[2]},
        voxel_types,
    };
    double x, y, z, grid_x, grid_y, grid_z;
    double delta[3] = { voxels.res_voxel[0] / 2.0, voxels.res_voxel[1] / 2.0, voxels.res_voxel[2] / 2.0};
    for (int i = 0; i < voxels.n_voxel[0]; i++)
    {
        grid_x = voxels.start_pt[0] + voxels.res_voxel[0]*i;
        x = grid_x + delta[0];
        voxels.centers.push_back(std::vector<std::vector<Eigen::Vector3d>>());
        voxels.voxel_pts.push_back(std::vector<std::vector<Eigen::Vector3d>>());
        for (int j = 0; j < voxels.n_voxel[1]; j++)
        {
            grid_y = voxels.start_pt[1] + voxels.res_voxel[1]*j;
            y = grid_y + delta[1];
            voxels.centers[i].push_back(std::vector<Eigen::Vector3d>());
            voxels.voxel_pts[i].push_back(std::vector<Eigen::Vector3d>());
            for (int k = 0; k < voxels.n_voxel[2]; k++) 
            {
                grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                z = grid_z + delta[2];
                Eigen::Vector3d center_p(x, y, z);
                Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                voxels.centers[i][j].push_back(center_p);
                voxels.voxel_pts[i][j].push_back(grid_pt);
                if (k == voxels.n_voxel[2]-1) 
                {
                    grid_z += voxels.res_voxel[2];
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[i][j].push_back(grid_pt);
                }
            }
            if (j == voxels.n_voxel[1]-1)
            {
                grid_y += voxels.res_voxel[1];
                voxels.voxel_pts[i].push_back(std::vector<Eigen::Vector3d>());
                for (int k = 0; k < voxels.n_voxel[2]; k++) 
                {
                    grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[i][j+1].push_back(grid_pt);
                    if (k == voxels.n_voxel[2]-1) 
                    {
                        grid_z += voxels.res_voxel[2];
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[i][j+1].push_back(grid_pt);
                    }
                }   
            }
        }
        if (i == voxels.n_voxel[0]-1)
        {
            grid_x += voxels.res_voxel[0];
            voxels.voxel_pts.push_back(std::vector<std::vector<Eigen::Vector3d>>());
            for (int j = 0; j < voxels.n_voxel[1]; j++)
            {
                grid_y = voxels.start_pt[1] + voxels.res_voxel[1]*j;
                voxels.voxel_pts[i+1].push_back(std::vector<Eigen::Vector3d>());
                for (int k = 0; k < voxels.n_voxel[2]; k++) 
                {
                    grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[i+1][j].push_back(grid_pt);
                    if (k == voxels.n_voxel[2]-1) 
                    {
                        grid_z += voxels.res_voxel[2];
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[i+1][j].push_back(grid_pt);
                    }
                }
                if (j == voxels.n_voxel[1]-1)
                {
                    grid_y += voxels.res_voxel[1];
                    voxels.voxel_pts[i+1].push_back(std::vector<Eigen::Vector3d>());
                    for (int k = 0; k < voxels.n_voxel[2]; k++) 
                    {
                        grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[i+1][j+1].push_back(grid_pt);
                        if (k == voxels.n_voxel[2]-1) 
                        {
                            grid_z += voxels.res_voxel[2];
                            Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                            voxels.voxel_pts[i+1][j+1].push_back(grid_pt);
                        }
                    }   
                }
            }
        }
    }
}


void identifyFeatureVoxels(BBVoxel& voxels)
{
    voxels.voxel_types.clear();
    

}



void generateCage(const Eigen::MatrixXd& mesh_vertices, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{   
    
    Eigen::MatrixXd pca_basic_matrix;
    Eigen::MatrixXd pca_based_mesh_vertices;
    Eigen::Vector3d barycenter;
    calculatePCA(mesh_vertices, pca_basic_matrix, pca_based_mesh_vertices, barycenter);

    computeBB(pca_based_mesh_vertices, cage_vertices, cage_faces, pca_basic_matrix, barycenter);

    // float sparse_factor = 0.3 // optinally user can define the sparse factor 
    BBVoxel voxels;
    voxelizeBB(mesh_vertices, cage_vertices, voxels);

    identifyFeatureVoxels(voxels);
}; 