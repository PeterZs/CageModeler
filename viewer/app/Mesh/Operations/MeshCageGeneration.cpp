#include <Mesh/Operations/MeshCageGeneration.h>

#include <cagegenerations/CoarseBoundingCagesByXian.h>

MeshCageGeneration::ExecutionResult MeshCageGeneration::Execute()
{
    LOG_DEBUG("Computing New Cage");

    // Need to distinguish several methods for cage generation here
    // At the moment only one method is implemented
    
    if (_params._cageGenerationType == CageGenerationMethod::CoarseBoundingCagesByXian)
    {
        generateCageCoarseBouding(_params._mesh._vertices, _params._mesh._faces, _params._cage._vertices, _params._cage._faces, _params._setting._sparseFactor, _params._setting._cageSmoothFactor);   
    }

    MeshCageGenerationResult result;
    result._cage = std::move(_params._cage);

    return result;
} 