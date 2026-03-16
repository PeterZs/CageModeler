# pragma once

#include <Mesh/Operations/MeshOperation.h>
#include <Mesh/Operations/MeshWeightsParams.h>

#include <UI/NewCageSetting.h>

#include <Eigen/Core>

enum class CageGenerationMethod : uint8_t
{
	CoarseBoundingCagesByXian
};

struct MeshCageGenerationParams
{
    MeshCageGenerationParams(
        const CageGenerationMethod generationType,
        EigenMesh mesh,
        EigenMesh cage,
        NewCageSetting setting)
        : _cageGenerationType(generationType)
        , _mesh(mesh)
        , _cage(cage)
        , _setting(setting)
    { }

    CageGenerationMethod _cageGenerationType;
    EigenMesh _mesh;
    EigenMesh _cage;
    NewCageSetting _setting;
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