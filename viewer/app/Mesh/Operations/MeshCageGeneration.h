# pragma once

#include <Mesh/Operations/MeshOperation.h>
#include <Mesh/Operations/MeshWeightsParams.h>

#include <Eigen/Core>

struct MeshCageGenerationParams
{
    MeshCageGenerationParams(
        EigenMesh mesh,
        EigenMesh cage)
        : _mesh(mesh)
        , _cage(cage)
    { }


    EigenMesh _mesh;
    EigenMesh _cage;
    
};

class MeshCageGeneration final : public MeshOperationTemplated<MeshCageGenerationParams, MeshCageGenerationResult>
{
    public: 
    using MeshOperationTemplated::MeshOperationTemplated;

	[[nodiscard]] std::string GetDescription() const override
	{
		return "Generating new cage";
	}

	ExecutionResult Execute();
};