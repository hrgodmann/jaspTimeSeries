import QtQuick
import QtQuick.Layouts
import JASP
import JASP.Controls

Section
{
	title: qsTr("Model")
	columns: 1
	property bool hasCovariates: false
	info: qsTr("The latent level follows a random walk on the log-odds scale. Set the priors for its initial probability, state-noise standard deviation and, when applicable, covariate coefficients.")

	Group
	{
		title: qsTr("Initial Probability")
		info: qsTr("Beta prior for the initial baseline probability. Alpha = beta = 1 gives a uniform prior. With covariates, this is the baseline when all centered covariates equal zero, not necessarily the probability at their first observed values.")

		DoubleField
		{
			name: "priorTheta1Shape1"
			label: qsTr("Alpha")
			defaultValue: 1
			min: 0
			inclusive: JASP.MaxOnly
		}

		DoubleField
		{
			name: "priorTheta1Shape2"
			label: qsTr("Beta")
			defaultValue: 1
			min: 0
			inclusive: JASP.MaxOnly
		}
	}

	Group
	{
		title: qsTr("State Noise")
		info: qsTr("Half-normal prior for the standard deviation of changes in the latent log-odds level. Larger prior scales allow larger changes between consecutive time points.")

		DoubleField
		{
			name: "priorSigmaSD"
			label: qsTr("Half-normal scale")
			defaultValue: 1
			min: 0
			inclusive: JASP.MaxOnly
		}
	}

	RadioButtonGroup
	{
		name: "priorBetaDistribution"
		title: qsTr("Covariate Coefficients")
		visible: hasCovariates
		info: qsTr("The selected zero-centered prior is applied independently to each standardized covariate coefficient, using the same settings for all coefficients. Each coefficient is on the log-odds scale and is constant over time.")

		GridLayout
		{
			columns: 2
			rowSpacing: jaspTheme.rowGroupSpacing
			columnSpacing: 0

			RadioButton
			{
				id: priorBetaNormal
				value: "normal"
				label: qsTr("Normal")
				info: qsTr("Normal prior centered at zero, with the specified standard deviation.")
			}
			DoubleField
			{
				name: "priorBetaNormalSD"
				label: qsTr("SD")
				enabled: priorBetaNormal.checked
				defaultValue: 1
				min: 0
				inclusive: JASP.MaxOnly
				info: qsTr("Positive standard deviation of the Normal coefficient prior.")
			}

			RadioButton
			{
				id: priorBetaStudentT
				value: "studentT"
				label: qsTr("Student t")
				checked: true
				info: qsTr("Student t prior centered at zero. This is the default coefficient prior, with three degrees of freedom by default.")
			}
			RowLayout
			{
				enabled: priorBetaStudentT.checked
				DoubleField
				{
					name: "priorBetaScale"
					label: qsTr("Scale")
					defaultValue: 0.5
					min: 0
					inclusive: JASP.MaxOnly
					info: qsTr("Positive scale of the Student t coefficient prior. For df greater than two, its standard deviation is scale multiplied by sqrt(df / (df - 2)). For df of two or less, the prior has no finite variance.")
				}
				DoubleField
				{
					name: "priorBetaDf"
					label: qsTr("df")
					defaultValue: 3
					min: 0
					inclusive: JASP.MaxOnly
					info: qsTr("Positive degrees of freedom of the Student t coefficient prior; fractional values are allowed. Smaller values give heavier tails. A value of one gives a Cauchy prior with the same scale.")
				}
			}

			RadioButton
			{
				id: priorBetaCauchy
				value: "cauchy"
				label: qsTr("Cauchy")
				info: qsTr("Cauchy prior centered at zero, with the specified scale. It has no finite mean or variance.")
			}
			DoubleField
			{
				name: "priorBetaCauchyScale"
				label: qsTr("Scale")
				enabled: priorBetaCauchy.checked
				defaultValue: 0.707
				min: 0
				inclusive: JASP.MaxOnly
				info: qsTr("Positive scale of the Cauchy coefficient prior. This is not a standard deviation.")
			}
		}
	}
}
