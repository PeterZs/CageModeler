#include <cagegenerations/CoarseBoundingCagesByXian.h>

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <limits>
#include <queue>
#include <map>

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

void minMax(const double v0, const double v1, const double v2, double& min, double& max)
{
    min = max = v0;
    if (v1 < min) min = v1;
    if (v1 > max) max = v1;
    if (v2 < min) min = v2;
    if (v2 > max) max = v2;
}


bool planeBoxOverlap(const Eigen::Vector3d normal, const double d, const Eigen::Vector3d maxbox)
{
    Eigen::Vector3d vmin = Eigen::Vector3d::Zero(1,3);
    Eigen::Vector3d vmax = Eigen::Vector3d::Zero(1,3);
    for (int i = 0; i < 3; i++) 
    {
        if (normal[i] > 0.0) 
        {
            vmin[i] = -maxbox[i];
            vmax[i] = maxbox[i];
        } 
        else 
        {
            vmin[i] = maxbox[i];
            vmax[i] = -maxbox[i];
        }
    }
    if (normal.dot(vmin)+d > 0.0) return false;
    if (normal.dot(vmax)+d >= 0.0) return true;
    return false;
}


bool axistest_x01(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p0 = a * v(0, 1) - b * v(0, 2);
    double p2 = a * v(2, 1) - b * v(2, 2);
    if (p0 < p2) 
    {
        min = p0;
        max = p2;
    }
    else 
    {
        min = p2;
        max = p0;
    }
    double rad = fa * halfbox[1] + fb * halfbox[2];
    if (min > rad || max < -rad) return false; 
    return true;
}

bool axistest_y02(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p0 = -a * v(0, 0) + b * v(0, 2);
    double p2 = -a * v(2, 0) + b * v(2, 2);
    if (p0 < p2) 
    {
        min = p0;
        max = p2;
    }
    else 
    {
        min = p2;
        max = p0;
    }
    double rad = fa * halfbox[0] + fb * halfbox[2];
    if (min > rad || max < -rad) return false; 
    return true;
}

bool axistest_z12(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p1 = a * v(1, 0) - b * v(1, 1);
    double p2 = a * v(2, 0) - b * v(2, 1);
    if (p2 < p1) 
    {
        min = p2;
        max = p1;
    }
    else 
    {
        min = p1;
        max = p2;
    }
    double rad = fa * halfbox[0] + fb * halfbox[1];
    if (min > rad || max < -rad) return false; 
    return true;
}

bool axistest_z0(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p0 = a * v(0, 0) - b * v(0, 1);
    double p1 = a * v(1, 0) - b * v(1, 1);
    if (p0 < p1) 
    {
        min = p0;
        max = p1;
    }
    else 
    {
        min = p1;
        max = p0;
    }
    double rad = fa * halfbox[0] + fb * halfbox[1];
    if (min > rad || max < -rad) return false; 
    return true;
}

bool axistest_x2(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p0 = a * v(0, 1) - b * v(0, 2);
    double p1 = a * v(1, 1) - b * v(1, 2);
    if (p0 < p1) 
    {
        min = p0;
        max = p1;
    }
    else 
    {
        min = p1;
        max = p0;
    }
    double rad = fa * halfbox[1] + fb * halfbox[2];
    if (min > rad || max < -rad) return false; 
    return true;
}

bool axistest_y1(const double a, const double b, const double fa, const double fb, const Eigen::Matrix3d& v, const Eigen::Vector3d& halfbox)
{
    double min, max;
    double p0 = -a * v(0, 0) + b * v(0, 2);
    double p1 = -a * v(1, 0) - b * v(1, 2);
    if (p0 < p1) 
    {
        min = p0;
        max = p1;
    }
    else 
    {
        min = p1;
        max = p0;
    }
    double rad = fa * halfbox[0] + fb * halfbox[2];
    if (min > rad || max < -rad) return false; 
    return true;
}


// Use the Fast 3D Triangle-Box Overlap Testing to identify feature voxels
void identifyFeatureVoxels(BBVoxel& voxels, const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXi& mesh_faces)
{
    voxels.voxel_types.clear();
    Eigen::Matrix3d triangle = Eigen::Matrix3d::Zero(3, 3);
    Eigen::Vector3d halfbox(voxels.res_voxel[0]/2, voxels.res_voxel[1]/2, voxels.res_voxel[2]/2);
    
    for (size_t k1 = 0; k1 < voxels.centers.size(); k1++) 
    {
        voxels.voxel_types.push_back(std::vector<std::vector<int8_t>>());
        for (size_t k2 = 0; k2 < voxels.centers[k1].size(); k2++) 
        {
            voxels.voxel_types[k1].push_back(std::vector<int8_t>());
            for (size_t k3 = 0; k3 < voxels.centers[k1][k2].size(); k3++)
            {
                voxels.voxel_types[k1][k2].push_back(-1);
                for (int i = 0; i < mesh_faces.rows(); i++)
                {
                    for (int j = 0; j < 3; j++) 
                    {
                        int idx = mesh_faces(i, j);
                        triangle.row(j) = mesh_vertices.row(idx) - voxels.centers[k1][k2][k3].transpose();
                    }
                    // first test: Test the AABB against the minimal AABB around the triangle.
                    double min, max;
                    minMax(triangle(0, 0), triangle(1, 0), triangle(2, 0), min, max);
                    if (min > halfbox[0] || max < -halfbox[0])
                    {
                        continue;
                    }
                    minMax(triangle(0, 1), triangle(1, 1), triangle(2, 1), min, max);
                    if (min > halfbox[1] || max < -halfbox[1])
                    {
                        continue;
                    }
                    minMax(triangle(0, 2), triangle(1, 2), triangle(2, 2), min, max);
                    if (min > halfbox[2] || max < -halfbox[2])
                    {
                        continue;
                    }

                    // second test: fast plane/AABB overlap test 
                    Eigen::Vector3d e0 = triangle.row(1) - triangle.row(0);
                    Eigen::Vector3d e1 = triangle.row(2) - triangle.row(1);
                    Eigen::Vector3d normal = e0.cross(e1);
                    double d = -normal.dot(triangle.row(0));

                    if (!planeBoxOverlap(normal, d, halfbox)) 
                    {
                        continue;
                    }

                    // third test: crossproduct(edge from triangle, {x,y,z}-directin) 
                    Eigen::Vector3d e2 = triangle.row(0) - triangle.row(2);

                    Eigen::Vector3d abs = e0.cwiseAbs();
                    if (!axistest_x01(e0[2], e0[1], abs[2], abs[1], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_y02(e0[2], e0[0], abs[2], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_z12(e0[1], e0[0], abs[1], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }

                    abs = e1.cwiseAbs();
                    if (!axistest_x01(e1[2], e1[1], abs[2], abs[1], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_y02(e1[2], e1[0], abs[2], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_z0(e1[1], e1[0], abs[1], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }

                    abs = e2.cwiseAbs();
                    if (!axistest_x2(e2[2], e2[1], abs[2], abs[1], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_y1(e2[2], e2[0], abs[2], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }
                    if (!axistest_z12(e2[1], e2[0], abs[1], abs[0], triangle, halfbox)) 
                    {
                        continue;
                    }

                    voxels.voxel_types[k1][k2][k3] = 0;
                    break;
                }
                // break;
            }
            // break;
        }
        // break;
    }
}


void renderFeatureVoxelHelper(BBVoxel voxels, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{
    int num_of_voxels = voxels.n_voxel[0] * voxels.n_voxel[1] * voxels.n_voxel[2];
    int num_of_feature_voxels = 0;
    int num_of_outer_voxels = 0;
    int num_of_inner_voxels = 0;
    for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++)  
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++)
            {
                if (voxels.voxel_types[k1][k2][k3] == 0)    
                {
                    num_of_feature_voxels++;           
                }
                if (voxels.voxel_types[k1][k2][k3] == 1)    
                {
                    num_of_outer_voxels++;           
                }   
                if (voxels.voxel_types[k1][k2][k3] == -1)    
                {
                    num_of_inner_voxels++;           
                }   
            }
        }
    }

    cage_vertices = Eigen::MatrixXd::Zero(num_of_feature_voxels * 8, 3);
    cage_faces = Eigen::MatrixXi::Zero(num_of_feature_voxels * 12, 3);
    int i = 0;
    for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++)  
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++)
            {
                if (voxels.voxel_types[k1][k2][k3] == 0)    
                {
                    cage_vertices.row(0 + i*8) = voxels.voxel_pts[k1][k2][k3];
                    cage_vertices.row(1 + i*8) = voxels.voxel_pts[k1][k2][k3+1];
                    cage_vertices.row(2 + i*8) = voxels.voxel_pts[k1][k2+1][k3];
                    cage_vertices.row(3 + i*8) = voxels.voxel_pts[k1][k2+1][k3+1];
                    cage_vertices.row(4 + i*8) = voxels.voxel_pts[k1+1][k2][k3];
                    cage_vertices.row(5 + i*8) = voxels.voxel_pts[k1+1][k2][k3+1];
                    cage_vertices.row(6 + i*8) = voxels.voxel_pts[k1+1][k2+1][k3];
                    cage_vertices.row(7 + i*8) = voxels.voxel_pts[k1+1][k2+1][k3+1];

                    cage_faces.row(0 + i*12) = Eigen::Vector3i(0 + i*8, 1 + i*8, 4 + i*8);
                    cage_faces.row(1 + i*12) = Eigen::Vector3i(1 + i*8, 5 + i*8, 4 + i*8);
                    cage_faces.row(2 + i*12) = Eigen::Vector3i(0 + i*8, 2 + i*8, 1 + i*8);
                    cage_faces.row(3 + i*12) = Eigen::Vector3i(1 + i*8, 2 + i*8, 3 + i*8);
                    cage_faces.row(4 + i*12) = Eigen::Vector3i(4 + i*8, 5 + i*8, 7 + i*8);
                    cage_faces.row(5 + i*12) = Eigen::Vector3i(4 + i*8, 7 + i*8, 6 + i*8);
                    cage_faces.row(6 + i*12) = Eigen::Vector3i(3 + i*8, 6 + i*8, 7 + i*8);
                    cage_faces.row(7 + i*12) = Eigen::Vector3i(2 + i*8, 6 + i*8, 3 + i*8);
                    cage_faces.row(8 + i*12)  = Eigen::Vector3i(1 + i*8, 7 + i*8, 5 + i*8);
                    cage_faces.row(9 + i*12)  = Eigen::Vector3i(1 + i*8, 3 + i*8, 7 + i*8);
                    cage_faces.row(10 + i*12) = Eigen::Vector3i(0 + i*8, 4 + i*8, 2 + i*8);
                    cage_faces.row(11 + i*12) = Eigen::Vector3i(2 + i*8, 4 + i*8, 6 + i*8);
                                    
                    i++;
                }   
            }
        }
    }


    std::cout << "num features: " << num_of_feature_voxels << std::endl;
    std::cout << "num outer: " << num_of_outer_voxels << std::endl;
    std::cout << "num inner: " << num_of_inner_voxels << std::endl;
    std::cout << "num voxels: " << num_of_voxels << std::endl;
    // std::cout << (double)num_of_feature_voxels / num_of_voxels << std::endl;

}


void identifyOuterVoxelsWithFillingAlgo(BBVoxel& voxels, const Eigen::Vector3i seed)
{
    if (voxels.voxel_types[seed[0]][seed[1]][seed[2]] == 0 || voxels.voxel_types[seed[0]][seed[1]][seed[2]] == 1)
    {
        return;
    }

    const std::array<Eigen::Vector3i, 6> neighbors = {{
        {  1,  0,  0 },
        { -1,  0,  0 },
        {  0,  1,  0 },
        {  0, -1,  0 },
        {  0,  0,  1 },
        {  0,  0, -1 }
    }};

    std::queue<Eigen::Vector3i> queue;
    queue.push(seed); 
    while (!queue.empty())
    {
        Eigen::Vector3i indices = queue.front();
        queue.pop();

        if (voxels.voxel_types[indices[0]][indices[1]][indices[2]] == 0 || voxels.voxel_types[indices[0]][indices[1]][indices[2]] == 1)
        {
            continue; // if already visited or feature voxel --> just continue
        }

        voxels.voxel_types[indices[0]][indices[1]][indices[2]] = 1; //mark as outer voxels

        for (const auto& n: neighbors)
        {   
            Eigen::Vector3i n_indices = indices + n;

            if (n_indices[0] < 0 || n_indices[1] < 0 || n_indices[2] < 0 || n_indices[0] >= voxels.n_voxel[0] || n_indices[1] >= voxels.n_voxel[1] || n_indices[2] >= voxels.n_voxel[2])
            {
                continue;
            }
            if (voxels.voxel_types[n_indices[0]][n_indices[1]][n_indices[2]] != 0 && voxels.voxel_types[n_indices[0]][n_indices[1]][n_indices[2]] != 1)
            {
                queue.push(n_indices);
            }
        }
    }
}



void extractSingleFace(BBVoxel& voxels, std::vector<Eigen::Vector3i>& outer_faces, std::vector<int>& outer_vertices, const int& v0, const int& v1, const int& v2, const int& v3)
{
    int i0, i1, i2, i3 = -1;
    for (int i = 0; i < outer_vertices.size(); i++)
    {   
        if (outer_vertices[i] == v0)
        {
            i0 = i;
        } 
        else if (outer_vertices[i] == v1)
        {
            i1 = i;
        }
        else if (outer_vertices[i] == v2)
        {
            i2 = i;
        }
        else if (outer_vertices[i] == v3)
        {
            i3 = i;
        }
    }
    if (i0 != -1) 
    {
        i0 = outer_vertices.size();
        outer_vertices.push_back(v0);
    }
    if (i1 != -1) 
    {
        i1 = outer_vertices.size();
        outer_vertices.push_back(v1);
    }
    if (i2 != -1) 
    {
        i2 = outer_vertices.size();
        outer_vertices.push_back(v2);
    }
    if (i3 != -1) 
    {
        i3 = outer_vertices.size();
        outer_vertices.push_back(v3);
    }
    outer_faces.push_back(Eigen::Vector3i(i0, i1, i2));
    outer_faces.push_back(Eigen::Vector3i(i1, i2, i3));
}



void extractOuterSurface(BBVoxel& voxels, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{
    std::vector<Eigen::Vector3i> outer_faces;
    std::vector<int> outer_vertices;
    int vertex_id = 0;
    int v0, v1, v2, v3;
    // int vertices_id = k1 + k2 * voxels.n_voxel[0] + k3 * voxels.n_voxel[0] * voxels.n_voxel[1];
    // std::cout << vertices_id << std::endl;
    // std::cout << voxels.voxel_pts.size() * voxels.voxel_pts[0].size() * voxels.voxel_pts[0][0].size() << std::endl;

    int num_voxel_dim_0 = voxels.voxel_pts.size();
    int num_voxel_dim_1 = voxels.voxel_pts[0].size();

    for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++)  
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++)
            {
                int t = k1 + k2 * (voxels.voxel_pts.size()) + k3 * (voxels.voxel_pts.size()) * (voxels.voxel_pts[0].size());
                // std::cout << t << std::endl;
                if (voxels.voxel_types[k1][k2][k3] == -1 || voxels.voxel_types[k1][k2][k3] == 1) 
                {
                    continue;
                }   
                // only extract if given voxel is a feature voxel
                v0 = k1 + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = k1 + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = k1 + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = k1 + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                if (k1 == 0) // edge --> extract the outer surface
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                } 
                else if (voxels.voxel_types[k1-1][k2][k3] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                v0 = k1 + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = (k1+1) + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = (k1+1) + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = k1 + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                if (k2 == 0)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                else if (voxels.voxel_types[k1][k2-1][k3] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                v0 = k1 + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = (k1+1) + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = (k1+1) + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = k1 + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                if (k3 == 0)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                else if (voxels.voxel_types[k1][k2-1][k3] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                v0 = (k1+1) + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = (k1+1) + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = (k1+1) + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = (k1+1) + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                if (k1 == voxels.n_voxel[0]-1) // edge --> extract the outer surface
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                } 
                else if (voxels.voxel_types[k1+1][k2][k3] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                v0 = k1 + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = (k1+1) + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = (k1+1) + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = k1 + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                if (k2 == voxels.n_voxel[1]-1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                else if (voxels.voxel_types[k1][k2+1][k3] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                v0 = k1 + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v1 = (k1+1) + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v2 = (k1+1) + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                v3 = k1 + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1;
                if (k3 == voxels.n_voxel[2]-1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
                else if (voxels.voxel_types[k1][k2][k3+1] == 1)
                {
                    extractSingleFace(voxels, outer_faces, outer_vertices, v0, v1, v2, v3);
                }
            }
        }
    }

    cage_vertices = Eigen::MatrixXd::Zero(outer_vertices.size(), 3);
    for (int i = 0; i < outer_vertices.size(); i++)
    {
        int k1 = outer_vertices[i] % voxels.voxel_pts.size();
        int k2 = (outer_vertices[i] / voxels.voxel_pts.size()) % voxels.voxel_pts[0].size();
        int k3 = outer_vertices[i] / (voxels.voxel_pts.size() * voxels.voxel_pts[0].size());
        // std::cout << k1 << "  " << k2 << "  " << k3 << std::endl;
        cage_vertices.row(i) = voxels.voxel_pts[k1][k2][k3];
    }
    cage_faces = Eigen::MatrixXi::Zero(outer_faces.size(), 3);
    for (int i = 0; i < outer_faces.size(); i++)
    {
        cage_faces.row(i) = outer_faces[i];
    }
}


void generateCage(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXi& mesh_faces, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{   
    Eigen::MatrixXd pca_basic_matrix;
    Eigen::MatrixXd pca_based_mesh_vertices;
    Eigen::Vector3d barycenter;
    calculatePCA(mesh_vertices, pca_basic_matrix, pca_based_mesh_vertices, barycenter);

    computeBB(pca_based_mesh_vertices, cage_vertices, cage_faces, pca_basic_matrix, barycenter);

    float sparse_factor = 0.8; // optinally user can define the sparse factor 
    BBVoxel voxels;
    voxelizeBB(mesh_vertices, cage_vertices, voxels, sparse_factor);

    identifyFeatureVoxels(voxels, mesh_vertices, mesh_faces);

    identifyOuterVoxelsWithFillingAlgo(voxels, Eigen::Vector3i::Zero(1,3));

    extractOuterSurface(voxels, cage_vertices, cage_faces);

    // renderFeatureVoxelHelper(voxels, cage_vertices, cage_faces);
}