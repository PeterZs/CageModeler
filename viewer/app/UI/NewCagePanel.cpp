#include <UI/NewCagePanel.h>
#include <UI/ProjectSettingsHelpers.h>
#include <Mesh/Operations/MeshOperationSystem.h>

NewCagePanel::NewCagePanel( const std::shared_ptr<ProjectModelData>& model,
                            const std::shared_ptr<MeshOperationSystem>& meshOperationSystem,
                            std::function<void ()> cancelDelegate,
                            std::function<void()> createDelegate)
	: _meshOperationSystem(meshOperationSystem)
	, _cancelDelegate(std::move(cancelDelegate))
	, _createDelegate(std::move(createDelegate))
	, _createButtonPressed(false)
	, _cancelButtonPressed(false)
{
    _modifiedProjectModel = *model;

    SetModel(model);
}

void NewCagePanel::Layout()
{
    if (_isModalVisible)
	{
		ImGui::OpenPopup("ProjectSettings");
	}

	const auto displaySize = ImGui::GetIO().DisplaySize;

	ImGui::SetNextWindowPos(ImVec2(displaySize.x / 2.0f, displaySize.y / 2.0f),
		ImGuiCond_Appearing,
		ImVec2(0.5f, 0.5f));

	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(20.0f, 20.0f));

	if (ImGui::BeginPopupModal("ProjectSettings", &_isModalVisible,
		ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoSavedSettings))
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
		if (ImGui::BeginTable("##SettingsTable", 2, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_SizingFixedSame))
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

				if (ImGui::BeginCombo("##Deformation", ProjecSettingsHelpers::DeformationMethodNames[_selectedDeformationTypeIndex], ImGuiComboFlags_HeightRegular))
				{
					for (auto i = 0; i < ProjecSettingsHelpers::DeformationMethodNames.size(); i++)
					{
						const auto isSelected = (_selectedDeformationTypeIndex == i);

						if (ImGui::Selectable(ProjecSettingsHelpers::DeformationMethodNames[i], isSelected))
						{
							_selectedDeformationTypeIndex = i;
						}

						if (isSelected)
						{
							ImGui::SetItemDefaultFocus();
						}
					}

                    _modifiedProjectModel._deformationType = static_cast<DeformationType>(_selectedDeformationTypeIndex);

					ImGui::EndCombo();
				}
			}

			ImGui::BeginDisabled(_modifiedProjectModel._deformationType != DeformationType::LBC);
			{
				ImGui::TableNextRow();
				{
					ImGui::TableSetColumnIndex(0);
					ImGui::TextEx("LBC Weighting Scheme");
					ImGui::SameLine();
					UIHelpers::HelpMarker("The weighting scheme for LBC, which controls the level of locality.");
					ImGui::SameLine();

					ImGui::TableSetColumnIndex(1);

					UIHelpers::SetRightAligned(125.0f);

					if (ImGui::BeginCombo("##WeightingScheme", ProjecSettingsHelpers::LBCWeightingSchemeNames[_selectedWeightingSchemeIndex], ImGuiComboFlags_HeightRegular))
					{
						for (auto i = 0; i < ProjecSettingsHelpers::LBCWeightingSchemeNames.size(); i++)
						{
							const auto isSelected = (_selectedWeightingSchemeIndex == i);

							if (ImGui::Selectable(ProjecSettingsHelpers::LBCWeightingSchemeNames[i], isSelected))
							{
								_selectedWeightingSchemeIndex = i;
							}

							if (isSelected)
							{
								ImGui::SetItemDefaultFocus();
							}
						}

						ImGui::EndCombo();
					}

					_modifiedProjectModel._LBCWeightingScheme = static_cast<LBC::DataSetup::WeightingScheme>(_selectedWeightingSchemeIndex);
				}
			}
			ImGui::EndDisabled();

			ImGui::Dummy(ImVec2(0.0f, 8.0f));

			ImGui::TableNextRow();
			{
				ImGui::TableSetColumnIndex(0);

				ImGui::PushFont(UIStyle::BoldFont);
				ImGui::SetWindowFontScale(1.05f);
				ImGui::TextEx("Somigliana Settings");
				ImGui::SetWindowFontScale(1.0f);
				ImGui::PopFont();
			}

			ImGui::Dummy(ImVec2(0.0f, 5.0f));

			ImGui::BeginDisabled(_modifiedProjectModel._deformationType != DeformationType::Somigliana);
			{
				ImGui::InputDouble("##SomigNuNewCage", &_model->_somigNu);
				ImGui::SameLine();
				UIHelpers::HelpMarker("The material parameter nu for somigliana deformer.");

				/*auto bulgingValue = _model->_somigBulging.load();
				if (ImGui::InputDouble("##SomigBulging", &bulgingValue)) {
					_model->_somigBulging = bulgingValue;
				}
				ImGui::SameLine();
				UIHelpers::HelpMarker("The bulging parameter gamma for somigliana deformer.");

				auto blendFactorValue = _model->_somigBlendFactor.load();
				if (ImGui::InputDouble("##SomigBlendFactor", &blendFactorValue))
				{
					_model->_somigBlendFactor = blendFactorValue;
				}
				ImGui::SameLine();
				UIHelpers::HelpMarker("The blending factor for somigliana deformer interpolating between local and global boundary conditions.");

				if (ImGui::BeginCombo("Bulging Type",
					ProjecSettingsHelpers::SomiglianaBulgingTypeNames[_selectedBulgingTypeIndex],
					ImGuiComboFlags_HeightRegular | ImGuiComboFlags_WidthFitPreview))
				{
					for (auto i = 0; i < ProjecSettingsHelpers::SomiglianaBulgingTypeNames.size(); i++)
					{
						const auto isSelected = (_selectedBulgingTypeIndex == i);

						if (ImGui::Selectable(ProjecSettingsHelpers::SomiglianaBulgingTypeNames[i], isSelected))
						{
							_selectedBulgingTypeIndex = i;
						}

						if (isSelected)
						{
							ImGui::SetItemDefaultFocus();
						}
					}

					ImGui::EndCombo();
				}

				_model->_somigBulgingType = static_cast<BulgingType>(_selectedBulgingTypeIndex);

				ImGui::SameLine();
				UIHelpers::HelpMarker("The bulging type for somigliana deformer.");*/
			}
			ImGui::EndDisabled();

			ImGui::Dummy(ImVec2(0.0f, 8.0f));

			ImGui::TableNextRow();
			{
				ImGui::TableSetColumnIndex(0);

				ImGui::PushFont(UIStyle::BoldFont);
				ImGui::SetWindowFontScale(1.05f);
				ImGui::TextEx("Cage Generation Settings");
				ImGui::SetWindowFontScale(1.0f);
				ImGui::PopFont();
			}

			ImGui::Dummy(ImVec2(0.0f, 5.0f));

            ImGui::TableNextRow();
            {
                ImGui::TableSetColumnIndex(0);
                ImGui::TextEx("Methods");
                ImGui::SameLine();

                ImGui::TableSetColumnIndex(1);

                UIHelpers::SetRightAligned(125.0f);

                if (ImGui::BeginCombo("##CageGeneration", ProjecSettingsHelpers::CageGenerationMethodNames[_selectedCageGenerationIndex], ImGuiComboFlags_HeightRegular))
                {
                    for (auto i = 0; i < ProjecSettingsHelpers::CageGenerationMethodNames.size(); i++)
                    {
                        const auto isSelected = (_selectedCageGenerationIndex == i);

                        if (ImGui::Selectable(ProjecSettingsHelpers::CageGenerationMethodNames[i], isSelected))
                        {
                            _selectedCageGenerationIndex = i;
                        }

                        if (isSelected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    
                    _modifiedProjectModel._cageGenerationType = static_cast<CageGenerationMethod>(_selectedCageGenerationIndex);

                    ImGui::EndCombo();
                }
            }

            ImGui::Dummy(ImVec2(0.0f, 5.0f));

            ImGui::TableNextRow();
            {
                ImGui::TableSetColumnIndex(0);
                ImGui::TextEx("Sparse Factor");
                ImGui::SameLine();
                UIHelpers::HelpMarker("Sparse Factor, scale cage size.");

                ImGui::TableSetColumnIndex(1);
                UIHelpers::SetRightAligned(100.0f);
                ImGui::InputFloat("##Factor", &_setting._sparseFactor, 0.01f, 0.0f, "%.2f", ImGuiInputTextFlags_NoHorizontalScroll);
            }

			ImGui::Dummy(ImVec2(0.0f, 5.0f));

            ImGui::TableNextRow();
            {
                ImGui::TableSetColumnIndex(0);
                ImGui::TextEx("Cage Smoothness Factor");
                ImGui::SameLine();
                UIHelpers::HelpMarker("Cage Smoothness Factor.");

                ImGui::TableSetColumnIndex(1);
                UIHelpers::SetRightAligned(100.0f);
                ImGui::InputFloat("##CageSmoothFactor", &_setting._cageSmoothFactor, 0.01f, 0.0f, "%.2f", ImGuiInputTextFlags_NoHorizontalScroll);
            }
		}
		ImGui::EndDisabled();

		ImGui::EndTable();

		ImGui::Dummy(ImVec2(0.0f, 15.0f));

		const auto buttonSize = ImGui::CalcItemSize(ImVec2(100.0f, 0.0f), 0.0f, 0.0f);
		ImGui::SetCursorPosX(descColumnWidth + settingsColumnWidth - buttonSize.x * 2.0f - ImGui::GetStyle().ItemSpacing.x + ImGui::GetStyle().WindowPadding.x);

		// Cancel the creation and dismiss the popup.
		if (ImGui::Button("Cancel", buttonSize))
		{
			// No need to clear anything on the modified data since this is going to destroy the object anyways.

			Dismiss();
		}

		ImGui::SameLine();

		const auto hasNoMesh = !_modifiedProjectModel._meshFilepath.has_value() || !_modifiedProjectModel._cageFilepath.has_value();
		const auto hasNoDeformedCage = !_modifiedProjectModel._cageFilepath.has_value();
		const auto hasNoEmbedding = (DeformationTypeHelpers::RequiresEmbedding(_modifiedProjectModel._deformationType) && !_modifiedProjectModel._embeddingFilepath.has_value());
		const auto hasNoChanges = (_modifiedProjectModel == *_model);

		// Create the new cage.
		if (ImGui::Button("Create", buttonSize))
        {
            // Update the value of the actual model pointer.
            *_model = _modifiedProjectModel;

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

	_selectedDeformationTypeIndex = static_cast<uint32_t>(_model->_deformationType);
	_selectedWeightingSchemeIndex = static_cast<uint32_t>(_model->_LBCWeightingScheme);
    _selectedCageGenerationIndex = static_cast<uint32_t>(_model->_cageGenerationType);

	//_selectedBulgingTypeIndex = static_cast<uint32_t>(_model->_somigBulgingType);

	_modifiedProjectModel = *_model;
}

std::shared_ptr<ProjectModelData> NewCagePanel::GetModel() const 
{ 
	return _model; 
}


NewCageSetting NewCagePanel::GetSetting() const
{
    return _setting;
}