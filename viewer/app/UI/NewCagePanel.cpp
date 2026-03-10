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
        const auto descColumnWidth = std::max(0.2f * displaySize.x, 150.0f);
		const auto settingsColumnWidth = std::max(0.3f * displaySize.x, 250.0f);
		const auto horizontalOffset = descColumnWidth + settingsColumnWidth;

		const auto meshOperationSystem = _meshOperationSystem.lock();
		if (meshOperationSystem == nullptr)
		{
			return;
		}

        const auto hasRunningOperation = (meshOperationSystem->GetCurrentOperation() != nullptr);

        ImGui::BeginDisabled(hasRunningOperation);
        if (ImGui::BeginTable("##NewCageTable", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_SizingFixedSame))
        {
            constexpr auto columnFlags = ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoResize | ImGuiTableColumnFlags_NoReorder | ImGuiTableColumnFlags_NoHeaderLabel;

			ImGui::TableSetupColumn("##Description", columnFlags, descColumnWidth);
			ImGui::TableSetupColumn("##Settings", columnFlags, settingsColumnWidth);
        
            ImGui::TableNextRow();
            {
                ImGui::TableSetColumnIndex(0);
				ImGui::TextEx("Coordinates");
				ImGui::SameLine();

				ImGui::TableSetColumnIndex(1);

				UIHelpers::SetRightAligned(125.0f);
            }
        }





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
	UIPanel::SetModel(model);

	// _selectedDeformationTypeIndex = static_cast<uint32_t>(_model->_deformationType);
	// _selectedWeightingSchemeIndex = static_cast<uint32_t>(_model->_LBCWeightingScheme);

	// //_selectedBulgingTypeIndex = static_cast<uint32_t>(_model->_somigBulgingType);

	// _modifiedProjectModel = *_model;
}
