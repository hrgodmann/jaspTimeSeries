import QtQuick
import JASP.Controls
import "qml_components" as TS

Form
{
	columns: 1
	info: qsTr("Estimate a changing probability from successes and trials using a binomial state space model. Each row represents one equally spaced time point. At least two rows are required; missing values are not supported. Estimates use all supplied observations and are not forecasts.")

	VariablesForm
	{
		AvailableVariablesList { name: "allVariablesList" }

		AssignedVariablesList
		{
			name: "successes"
			title: qsTr("Successes")
			allowedColumns: ["scale"]
			singleVariable: true
			info: qsTr("Number of successes or events at each time point. Values must be nonnegative whole numbers and cannot exceed the corresponding number of trials. For binary observations, use 0 or 1 and a constant of one trial.")
		}

		DropDown
		{
			id: trialsMode
			name: "trialsMode"
			label: qsTr("Number of trials")
			startValue: "variable"
			values: [
				{ label: qsTr("Variable"), value: "variable" },
				{ label: qsTr("Constant"), value: "constant" }
			]
			info: qsTr("Use a variable when the number of trials differs between time points, or a constant when it is the same for every row. Trials are the observed opportunities for a success, not a population size.")
		}

		AssignedVariablesList
		{
			name: "trials"
			title: qsTr("Trials")
			allowedColumns: ["scale"]
			singleVariable: true
			visible: trialsMode.currentValue === "variable"
			info: qsTr("Number of trials at each time point, as nonnegative whole numbers. A row with zero trials and zero successes represents an unobserved time step. At least one row must have positive trials.")
		}

		IntegerField
		{
			name: "trialsConstant"
			label: qsTr("Trials per time point")
			defaultValue: 1
			min: 1
			max: 2147483647
			visible: trialsMode.currentValue === "constant"
			info: qsTr("The same number of trials is used for every row. Set this to one for binary observations.")
		}

		AssignedVariablesList
		{
			name: "time"
			title: qsTr("Time (optional)")
			allowedColumns: ["scale"]
			singleVariable: true
			info: qsTr("Optional numeric time variable. Rows are sorted by time; values must be unique and equally spaced. If left empty, row order defines equally spaced time steps. Keep unobserved steps as rows with zero successes and zero trials.")
		}

		// Keep list-valued covariate options. Both
		// singleVariable and maxRows: 1 would switch JASP to scalar options.
		AssignedVariablesList
		{
			id: covariates
			name: "covariates"
			title: qsTr("Covariates")
			allowedColumns: ["scale"]
			info: qsTr("Optional continuous covariates with additive effects that are constant over time. Leave empty for a model without covariates. Each covariate must vary across time and is separately centered and divided by its sample standard deviation. Covariates must not be linearly dependent after centering. Coefficients describe changes in log-odds, holding the other covariates and latent level fixed.")
		}
	}

	TS.BinomialStateSpacePlots
	{
		hasCovariates: covariates.count > 0
	}

	TS.BinomialStateSpaceModel
	{
		hasCovariates: covariates.count > 0
	}

	TS.BinomialStateSpaceDiagnostics
	{
		hasCovariates: covariates.count > 0
	}

	TS.BinomialStateSpaceAdvanced {}
}
