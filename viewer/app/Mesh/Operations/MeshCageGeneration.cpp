#include <Mesh/Operations/MeshCageGeneration.h>

#include <cagegenerations/CoarseBoundingCagesByXian.h>

MeshCageGeneration::ExecutionResult MeshCageGeneration::Execute()
{
    LOG_DEBUG("Computing New Cage");

    // Need to distinguish several methods for cage generation here
    // At the moment only one method is implemented

    generateCage(_params._mesh._vertices, _params._mesh._faces, _params._cage._vertices, _params._cage._faces);

    MeshCageGenerationResult result;
    result._cage_faces = std::move(_params._cage._faces);
    result._cage_vertices = std::move(_params._cage._vertices);

    return result;
}