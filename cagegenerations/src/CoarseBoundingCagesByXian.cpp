#include <cagegenerations/CoarseBoundingCagesByXian.h>

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <igl/readPLY.h>
#include <igl/cotmatrix.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/Hit.h>

#include <limits>
#include <set>
#include <queue>
#include <map>

static const double EPSILON_BB = 0.8;
static const float PREDEFINED_SPARSE_FACTOR = 0.1;


/**
 * Performs Principal Component Analysis (PCA) on mesh vertices.
 * Centers the data, computes the covariance matrix, and projects the vertices 
 * into the local coordinate system defined by the eigenvectors.
 */
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


/**
 * Computes an Oriented Bounding Box (OBB) for the mesh.
 * The function determines the axis-aligned min/max bounds in the PCA-transformed 
 * space, applies a small epsilon padding, generates a box mesh (8 vertices, 12 faces), 
 * and transforms the resulting cage back into the original world space.
 */
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
}


/**
 * Voxelizes the Bounding Box of the mesh.
 * This function calculates the voxel grid dimensions and resolution based on a 
 * sparseness factor and the mesh density. It initializes the 3D grid, computing 
 * both the center points of each voxel and the corner grid points.
 */
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
    std::vector<Eigen::Vector3d> centers;
    std::vector<Eigen::Vector3d> voxel_pts;
    std::vector<int8_t> voxel_types;
    voxels = BBVoxel{
        min_coordinates,
        centers,
        voxel_pts,
        {num_of_voxel_each_dim[0], num_of_voxel_each_dim[1], num_of_voxel_each_dim[2]},
        {res_of_voxel_each_dim[0], res_of_voxel_each_dim[1], res_of_voxel_each_dim[2]},
        voxel_types,
    };
    voxels.centers.resize(voxels.n_voxel[0] * voxels.n_voxel[1] * voxels.n_voxel[2]);
    voxels.voxel_pts.resize((voxels.n_voxel[0]+1) * (voxels.n_voxel[1]+1) * (voxels.n_voxel[2]+1));

    double x, y, z, grid_x, grid_y, grid_z;
    double delta[3] = { voxels.res_voxel[0] / 2.0, voxels.res_voxel[1] / 2.0, voxels.res_voxel[2] / 2.0};
    for (int i = 0; i < voxels.n_voxel[0]; i++)
    {
        grid_x = voxels.start_pt[0] + voxels.res_voxel[0]*i;
        x = grid_x + delta[0];
        for (int j = 0; j < voxels.n_voxel[1]; j++)
        {
            grid_y = voxels.start_pt[1] + voxels.res_voxel[1]*j;
            y = grid_y + delta[1];
            for (int k = 0; k < voxels.n_voxel[2]; k++) 
            {
                grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                z = grid_z + delta[2];
                Eigen::Vector3d center_p(x, y, z);
                Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                voxels.centers[voxels.getVoxelId(i, j, k)] = center_p;
                voxels.voxel_pts[voxels.getVoxelPtsId(i, j, k)] = grid_pt;
                if (k == voxels.n_voxel[2]-1) 
                {
                    grid_z += voxels.res_voxel[2];
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[voxels.getVoxelPtsId(i, j, k+1)] = grid_pt;
                }
            }
            if (j == voxels.n_voxel[1]-1)
            {
                grid_y += voxels.res_voxel[1];
                for (int k = 0; k < voxels.n_voxel[2]; k++) 
                {
                    grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[voxels.getVoxelPtsId(i, j+1, k)] = grid_pt;
                    if (k == voxels.n_voxel[2]-1) 
                    {
                        grid_z += voxels.res_voxel[2];
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[voxels.getVoxelPtsId(i, j+1, k+1)] = grid_pt;
                    }
                }   
            }
        }
        if (i == voxels.n_voxel[0]-1)
        {
            grid_x += voxels.res_voxel[0];
            for (int j = 0; j < voxels.n_voxel[1]; j++)
            {
                grid_y = voxels.start_pt[1] + voxels.res_voxel[1]*j;
                for (int k = 0; k < voxels.n_voxel[2]; k++) 
                {
                    grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                    Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                    voxels.voxel_pts[voxels.getVoxelPtsId(i+1, j, k)] = grid_pt;
                    if (k == voxels.n_voxel[2]-1) 
                    {
                        grid_z += voxels.res_voxel[2];
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[voxels.getVoxelPtsId(i+1, j, k+1)] = grid_pt;
                    }
                }
                if (j == voxels.n_voxel[1]-1)
                {
                    grid_y += voxels.res_voxel[1];
                    for (int k = 0; k < voxels.n_voxel[2]; k++) 
                    {
                        grid_z = voxels.start_pt[2] + voxels.res_voxel[2]*k;
                        Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                        voxels.voxel_pts[voxels.getVoxelPtsId(i+1, j+1, k)] = grid_pt;
                        if (k == voxels.n_voxel[2]-1) 
                        {
                            grid_z += voxels.res_voxel[2];
                            Eigen::Vector3d grid_pt(grid_x, grid_y, grid_z);
                            voxels.voxel_pts[voxels.getVoxelPtsId(i+1, j+1, k+1)] = grid_pt;
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

    for (size_t k1 = 0; k1 < voxels.n_voxel[0]; k1++)
    {
        for (size_t k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (size_t k3 = 0; k3 < voxels.n_voxel[2]; k3++)
            {
                voxels.voxel_types.push_back(-1);
                for (int i = 0; i < mesh_faces.rows(); i++)
                {
                    for (int j = 0; j < 3; j++) 
                    {
                        int idx = mesh_faces(i, j);
                        triangle.row(j) = mesh_vertices.row(idx) - voxels.centers[voxels.getVoxelId(k1, k2, k3)].transpose();
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

                    voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] = 0;
                    break;
                }
            }
        }
    }
}


void identifyOuterVoxels(BBVoxel& voxels)
{
    for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++)
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++)
            {
                if (k3 == 0 || k2 == 0 || k1 == 0 || k3 == voxels.n_voxel[2]-1 || k2 == voxels.n_voxel[1]-1 || k1 == voxels.n_voxel[0]-1)
                {
                    Eigen::Vector3i seed(k1, k2, k3);
                    identifyOuterVoxelsWithFillingAlgo(voxels, seed);
                }
            }
        }
    }
}


void identifyOuterVoxelsWithFillingAlgo(BBVoxel& voxels, const Eigen::Vector3i seed)
{
    if (voxels.voxel_types[voxels.getVoxelId(seed[0], seed[1], seed[2])] == 0 || voxels.voxel_types[voxels.getVoxelId(seed[0], seed[1], seed[2])] == 1)
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

        if (voxels.voxel_types[voxels.getVoxelId(indices[0], indices[1], indices[2])] == 0 || voxels.voxel_types[voxels.getVoxelId(indices[0], indices[1], indices[2])] == 1)
        {
            continue; // if already visited or feature voxel --> just continue
        }

        voxels.voxel_types[voxels.getVoxelId(indices[0], indices[1], indices[2])] = 1;

        for (const auto& n: neighbors)
        {   
            Eigen::Vector3i n_indices = indices + n;

            if (n_indices[0] < 0 || n_indices[1] < 0 || n_indices[2] < 0 || n_indices[0] >= voxels.n_voxel[0] || n_indices[1] >= voxels.n_voxel[1] || n_indices[2] >= voxels.n_voxel[2])
            {
                continue;
            }
            if (voxels.voxel_types[voxels.getVoxelId(n_indices[0], n_indices[1], n_indices[2])] != 0 && voxels.voxel_types[voxels.getVoxelId(n_indices[0], n_indices[1], n_indices[2])] != 1)
            {
                queue.push(n_indices);
            }
        }
    }
}


/**
 * Extracts a single face (quad) from the voxel grid and decomposes it into two triangles.
 * It checks if the provided vertex indices already exist in the 'outer_vertices' list
 * to avoid duplicates, updates the list if necessary, and stores the resulting 
 * triangle faces with consistent winding order.
 */
void extractSingleFace(std::vector<Eigen::Vector3i>& outer_faces, std::vector<int>& outer_vertices, const int& v0, const int& v1, const int& v2, const int& v3)
{
    int i0 = -1;
    int i1 = -1;
    int i2 = -1;
    int i3 = -1;
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
    if (i0 == -1) 
    {
        i0 = outer_vertices.size();
        outer_vertices.push_back(v0);
    }
    if (i1 == -1) 
    {
        i1 = outer_vertices.size();
        outer_vertices.push_back(v1);
    }
    if (i2 == -1) 
    {
        i2 = outer_vertices.size();
        outer_vertices.push_back(v2);
    }
    if (i3 == -1) 
    {
        i3 = outer_vertices.size();
        outer_vertices.push_back(v3);
    }
    // outer_faces.push_back(Eigen::Vector3i(i0, i1, i2));
    // outer_faces.push_back(Eigen::Vector3i(i2, i1, i3));
    outer_faces.push_back(Eigen::Vector3i(i0, i2, i1));
    outer_faces.push_back(Eigen::Vector3i(i1, i2, i3));
}


/**
 * Identifies and repairs non-manifold configurations within the voxel grid.
 * The function first resolves non-manifold edges by converting adjacent outer voxels 
 * into feature voxels (voxel attaching). In the second phase, it detects 
 * non-manifold vertices (where multiple regions meet at a single point) and 
 * marks them for later vertex splitting to ensure a topologically valid mesh.
 */
void checkForNonManifoldFeatureVoxels(BBVoxel& voxels)
{
    // for non manifold edge, simply perform voxel attaching
    // for non manifold vertex, more work has to be done.    
    for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++) 
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++) 
            {
                if (voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == 1 || voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == -1)
                {
                    continue;
                }
                // for non manifold egde
                // mark an outer voxel as feature to solve this
                int i1 = k1 - 1;
                int i2 = k2 - 1;
                int i3 = k3 - 1;
                int j1 = k1 + 1;
                int j2 = k2 + 1;
                int j3 = k3 + 1;
                bool manifold_mesh = false;
                if (i1 >= 0 && i2 >= 0 && voxels.voxel_types[voxels.getVoxelId(i1, i2, k3)] == 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] = 0;
                    }
                }
                if (i1 >= 0 && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] = 0;
                    }
                }
                if (i2 >= 0 && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] = 0;
                    }
                }
                if (j1 < voxels.n_voxel[0] && i2 >= 0 && voxels.voxel_types[voxels.getVoxelId(j1, i2, k3)] == 0 && voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] = 0;
                    }
                }
                if (j1 < voxels.n_voxel[0] && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(j1, k2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] = 0;
                    }
                }
                if (j2 < voxels.n_voxel[1] && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(k1, j2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(k1, j2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    manifold_mesh = true;
                    if (voxels.voxel_types[voxels.getVoxelId(k1, j2, k3)] == 1)
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, j2, k3)] = 0;
                    } 
                    else 
                    {
                        voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] = 0;
                    }
                }
                if (manifold_mesh)
                {
                    k1 = 0;
                    k2 = 0;
                    k3 = 0;
                }
            }
        }
    }

    // after fixing non manifold edges, look for non manifold vertices
    // in that case, mark them (vertices) as non manifold
    // vertex splitting needs to be executed but no in this method
    int num_voxel_ptsdim_0 = voxels.n_voxel[0] + 1;
    int num_voxel_ptsdim_1 = voxels.n_voxel[1] + 1;
    voxels.non_manifold_vertices.clear();
    voxels.splitted_vertices.clear();
    for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++) 
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++) 
            {
                if (voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == 1 || voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == -1)
                {
                    continue;
                }
                // for non manifold egde
                // mark an outer voxel as feature to solve this
                int i1 = k1 - 1;
                int i2 = k2 - 1;
                int i3 = k3 - 1;
                int j1 = k1 + 1;
                int j2 = k2 + 1;
                int j3 = k3 + 1;
                if (i1 >= 0 && i2 >= 0 && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(i1, i2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    int non_manifold_vertex = k1 + k2 * num_voxel_ptsdim_0 + k3 * num_voxel_ptsdim_0 * num_voxel_ptsdim_1;
                    voxels.non_manifold_vertices.push_back(non_manifold_vertex);
                    voxels.splitted_vertices.push_back(voxels.voxel_pts[voxels.getVoxelPtsId(k1, k2, k3)]);
                }
                if (j1 < voxels.n_voxel[0] && i2 >= 0 && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(j1, i2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(j1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0) 
                {
                    int non_manifold_vertex = (k1+1) + k2 * num_voxel_ptsdim_0 + k3 * num_voxel_ptsdim_0 * num_voxel_ptsdim_1;
                    voxels.non_manifold_vertices.push_back(non_manifold_vertex);
                    voxels.splitted_vertices.push_back(voxels.voxel_pts[voxels.getVoxelPtsId(k1+1, k2, k3)]);
                }
                if (i1 >= 0 && j2 < voxels.n_voxel[1] && i3 >= 0 && voxels.voxel_types[voxels.getVoxelId(i1, j2, i3)] == 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, j2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, i3)] != 0)
                {
                    int non_manifold_vertex = k1 + (k2+1) * num_voxel_ptsdim_0 + k3 * num_voxel_ptsdim_0 * num_voxel_ptsdim_1;
                    voxels.non_manifold_vertices.push_back(non_manifold_vertex);
                    voxels.splitted_vertices.push_back(voxels.voxel_pts[voxels.getVoxelPtsId(k1, k2+1, k3)]);
                }
                if (i1 >= 0 && i2 >= 0 && j3 >= 0 && voxels.voxel_types[voxels.getVoxelId(i1, i2, j3)] == 0 && voxels.voxel_types[voxels.getVoxelId(i1, k2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, i2, k3)] != 0 && voxels.voxel_types[voxels.getVoxelId(k1, k2, j3)] != 0)
                {
                    int non_manifold_vertex = k1 + k2 * num_voxel_ptsdim_0 + (k3+1) * num_voxel_ptsdim_0 * num_voxel_ptsdim_1;
                    voxels.non_manifold_vertices.push_back(non_manifold_vertex);
                    voxels.splitted_vertices.push_back(voxels.voxel_pts[voxels.getVoxelPtsId(k1, k2, k3+1)]);
                }
            }
        }
    }
}


void extractOuterSurface(BBVoxel& voxels, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces)
{
    // first check for non manifoldness in CBC.
    checkForNonManifoldFeatureVoxels(voxels);

    std::vector<Eigen::Vector3i> outer_faces;
    std::vector<int> outer_vertices;
    int vertex_id = 0;

    int num_voxel_dim_0 = voxels.n_voxel[0] + 1;
    int num_voxel_dim_1 = voxels.n_voxel[1] + 1;

    int num_of_all_voxel_pts = (voxels.n_voxel[0]+1) * (voxels.n_voxel[1]+1) * (voxels.n_voxel[2]+1);

    std::vector<bool> non_manifold_vertex_is_splitted(voxels.non_manifold_vertices.size(), false);  

    for (int k3 = 0; k3 < voxels.n_voxel[2]; k3++)  
    {
        for (int k2 = 0; k2 < voxels.n_voxel[1]; k2++)
        {
            for (int k1 = 0; k1 < voxels.n_voxel[0]; k1++)
            {               
                if (voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == -1 || voxels.voxel_types[voxels.getVoxelId(k1, k2, k3)] == 1)
                {
                    continue;
                }
                // only extract if given voxel is a feature voxel
                int voxel_vertices[8] = {
                    k1 + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1,
                    k1 + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1,
                    k1 + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1,
                    k1 + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1,
                    (k1+1) + k2 * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1,
                    (k1+1) + k2 * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1,
                    (k1+1) + (k2+1) * num_voxel_dim_0 + k3 * num_voxel_dim_0 * num_voxel_dim_1,
                    (k1+1) + (k2+1) * num_voxel_dim_0 + (k3+1) * num_voxel_dim_0 * num_voxel_dim_1
                };
                for (int i = 0; i < 8; i++)
                {
                    for (int j = 0; j < voxels.non_manifold_vertices.size(); j++) 
                    {
                        if (voxel_vertices[i] == voxels.non_manifold_vertices[j] && !non_manifold_vertex_is_splitted[j])
                        {
                            voxel_vertices[i] = num_of_all_voxel_pts + j;
                            non_manifold_vertex_is_splitted[j] = true;
                            continue;
                        }
                    }
                }

                if (k1 == 0) // edge --> extract the outer surface
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[0], voxel_vertices[2], voxel_vertices[1], voxel_vertices[3]);
                } 
                else if (voxels.voxel_types[voxels.getVoxelId(k1-1, k2, k3)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[0], voxel_vertices[2], voxel_vertices[1], voxel_vertices[3]);
                }
                if (k2 == 0)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[4], voxel_vertices[0], voxel_vertices[5], voxel_vertices[1]);
                }
                else if (voxels.voxel_types[voxels.getVoxelId(k1, k2-1, k3)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[4], voxel_vertices[0], voxel_vertices[5], voxel_vertices[1]);
                }
                if (k3 == 0)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[0], voxel_vertices[4], voxel_vertices[2], voxel_vertices[6]);
                }
                else if (voxels.voxel_types[voxels.getVoxelId(k1, k2, k3-1)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[0], voxel_vertices[4], voxel_vertices[2], voxel_vertices[6]);
                }
                if (k1 == voxels.n_voxel[0]-1) // edge --> extract the outer surface
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[4], voxel_vertices[5], voxel_vertices[6], voxel_vertices[7]);
                } 
                else if (voxels.voxel_types[voxels.getVoxelId(k1+1, k2, k3)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[4], voxel_vertices[5], voxel_vertices[6], voxel_vertices[7]);
                }
                if (k2 == voxels.n_voxel[1]-1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[2], voxel_vertices[6], voxel_vertices[3], voxel_vertices[7]);
                }
                else if (voxels.voxel_types[voxels.getVoxelId(k1, k2+1, k3)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[2], voxel_vertices[6], voxel_vertices[3], voxel_vertices[7]);
                }
                if (k3 == voxels.n_voxel[2]-1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[5], voxel_vertices[1], voxel_vertices[7], voxel_vertices[3]);
                }
                else if (voxels.voxel_types[voxels.getVoxelId(k1, k2, k3+1)] == 1)
                {
                    extractSingleFace(outer_faces, outer_vertices, voxel_vertices[5], voxel_vertices[1], voxel_vertices[7], voxel_vertices[3]);
                }
            }
        }
    }    

    cage_vertices = Eigen::MatrixXd::Zero(outer_vertices.size(), 3);
    for (size_t i = 0; i < outer_vertices.size(); i++)
    {
        if (outer_vertices[i] >= num_of_all_voxel_pts)
        {
            cage_vertices.row(i) = voxels.splitted_vertices[outer_vertices[i] - num_of_all_voxel_pts];
        }
        else 
        {
            int k1 = outer_vertices[i] % (voxels.n_voxel[0]+1);
            int k2 = (outer_vertices[i] / (voxels.n_voxel[0]+1)) % (voxels.n_voxel[1]+1);
            int k3 = outer_vertices[i] / ((voxels.n_voxel[0]+1) * (voxels.n_voxel[1]+1));
            cage_vertices.row(i) = voxels.voxel_pts[voxels.getVoxelPtsId(k1, k2, k3)];
        }
    }
    cage_faces = Eigen::MatrixXi::Zero(outer_faces.size(), 3);
    for (size_t i = 0; i < outer_faces.size(); i++)
    {
        cage_faces.row(i) = outer_faces[i];
    }
}



//////////////////////////////////7



void smoothCage(Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, float lambda_smooth, int it_smooth, Eigen::MatrixXd& pca_based_mesh_vertices, const Eigen::MatrixXi& mesh_faces)
{
    Eigen::SparseMatrix<double> laplace_beltrami_matrix_;
    igl::cotmatrix(cage_vertices,cage_faces,laplace_beltrami_matrix_);
    Eigen::MatrixXd cage_vertices_prev = cage_vertices;
    Eigen::MatrixXd cage_vertices_smooth = cage_vertices;
    Eigen::SparseMatrix<double> S(cage_vertices.rows(), cage_vertices.rows());
    S.setIdentity();
    Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solverImpl;
    S = S - lambda_smooth * laplace_beltrami_matrix_;
    solverImpl.compute(S);

    Eigen::Vector3d s, dir, dir_norm;
    std::vector<igl::Hit<double>> hits; 
    const double delta_T = 0.5;
    for (int i = 0; i < 2; i++){
        cage_vertices_smooth =  solverImpl.solve(cage_vertices_prev);
        // Test for intersection with the mesh

        for (int j = 0; j < cage_vertices.rows(); j++) {
            hits.clear();
            s = cage_vertices_prev.row(j);
            dir = cage_vertices_smooth.row(j) - cage_vertices_prev.row(j);
            dir_norm = dir.normalized();
            if (igl::ray_mesh_intersect(s,dir_norm,pca_based_mesh_vertices,mesh_faces,hits) || igl::ray_mesh_intersect(s,-dir_norm,pca_based_mesh_vertices,mesh_faces,hits)){
                if (hits[0].t < dir.norm()) {
                    cage_vertices_smooth.row(j) = Eigen::Vector3d(cage_vertices_prev.row(j)) + delta_T * hits[0].t * dir_norm;
                }
            }
        }
        cage_vertices_prev = cage_vertices_smooth;
    }

    cage_vertices = cage_vertices_smooth;   
}


void generateCageCoarseBouding(const Eigen::MatrixXd& mesh_vertices, const Eigen::MatrixXi& mesh_faces, Eigen::MatrixXd& cage_vertices, Eigen::MatrixXi& cage_faces, const float sparse_factor, const float cage_smooth_factor)
{   
    Eigen::MatrixXd pca_basic_matrix;
    Eigen::MatrixXd pca_based_mesh_vertices;
    Eigen::Vector3d barycenter;
    calculatePCA(mesh_vertices, pca_basic_matrix, pca_based_mesh_vertices, barycenter);

    computeBB(pca_based_mesh_vertices, cage_vertices, cage_faces, pca_basic_matrix, barycenter);
    
    BBVoxel voxels;
    voxelizeBB(mesh_vertices, cage_vertices, voxels, sparse_factor);

    identifyFeatureVoxels(voxels, mesh_vertices, mesh_faces);

    identifyOuterVoxels(voxels);

    extractOuterSurface(voxels, cage_vertices, cage_faces);

    float lambda_smooth = .0;
    if (cage_smooth_factor > 0.0)
    {
        lambda_smooth = cage_smooth_factor;
    }
    int num_it = 10;
    smoothCage(cage_vertices, cage_faces, lambda_smooth, num_it, pca_based_mesh_vertices, mesh_faces);
}