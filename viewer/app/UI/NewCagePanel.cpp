#include <UI/NewCagePanel.h>
#include <UI/ProjectSettingsHelpers.h>
#include <Mesh/Operations/MeshOperationSystem.h>

NewCagePanel::NewCagePanel( const std::shared_ptr<MeshOperationSystem>& meshOperationSystem,
                            std::function<void ()> cancelDelegate,
                            std::function<void()> createDelegate)
	: _meshOperationSystem(meshOperationSystem)
	, _cancelDelegate(std::move(cancelDelegate))
	, _createDelegate(std::move(createDelegate))
	, _createButtonPressed(false)
	, _cancelButtonPressed(false)
{

}

void NewCagePanel::Layout()
{
    if (!_isModalVisible)
        return;

    ImGui::OpenPopup("New Cage");

    const auto displaySize = ImGui::GetIO().DisplaySize;
    ImGui::SetNextWindowPos(ImVec2(displaySize.x / 2.0f, displaySize.y / 2.0f),
                             ImGuiCond_Appearing,
                             ImVec2(0.5f, 0.5f));

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(20.0f, 20.0f));

    if (ImGui::BeginPopupModal("New Cage", &_isModalVisible,
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove))
    {
        ImGui::Text("Create a new cage");

        ImGui::Dummy(ImVec2(0.0f, 10.0f)); 

        const auto buttonSize = ImVec2(100.0f, 0.0f);

        // Cancel Button
        if (ImGui::Button("Cancel", buttonSize))
        {
            Dismiss(); 
        }

        ImGui::SameLine();

        // Create Button
        if (ImGui::Button("Create", buttonSize))
        {
            _createButtonPressed = true;
        }

        ImGui::EndPopup();
    }

    ImGui::PopStyleVar();

    if (_createButtonPressed)
    {
        _createDelegate();
        _createButtonPressed = false;
    }
}

void NewCagePanel::Present()
{
	_isModalVisible = true;
}

void NewCagePanel::Dismiss()
{
	if (_isModalVisible)
	{
		ImGui::CloseCurrentPopup();

		_isModalVisible = false;
	}
}

void NewCagePanel::SetModel(const std::shared_ptr<ProjectModelData>& model)
{
	// UIPanel::SetModel(model);

	// _selectedDeformationTypeIndex = static_cast<uint32_t>(_model->_deformationType);
	// _selectedWeightingSchemeIndex = static_cast<uint32_t>(_model->_LBCWeightingScheme);

	// //_selectedBulgingTypeIndex = static_cast<uint32_t>(_model->_somigBulgingType);

	// _modifiedProjectModel = *_model;
}
