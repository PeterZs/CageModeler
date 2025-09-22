#include <Mesh/Operations/MeshCageGeneration.h>

#include <cagegenerations/CoarseBoundingCagesByXian.h>

MeshCageGeneration::ExecutionResult MeshCageGeneration::Execute()
{
    LOG_DEBUG("Computing New Cage");

    generateCage();

    return MeshCageGenerationResult {
        
    };
}