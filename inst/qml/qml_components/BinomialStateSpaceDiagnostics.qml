import QtQuick
import JASP.Controls

Section
{
	title: qsTr("MCMC Diagnostics")
	columns: 1
	property bool hasCovariates: false

	CheckBox
	{
		name: "showMcmcSummary"
		label: qsTr("Overview table")
		info: qsTr("Show posterior means, standard deviations, effective sample sizes and R-hat for the selected parameters. Review convergence before interpreting the estimates.")

		CheckBox
		{
			name: "showTheta"
			label: qsTr("Probabilities at all time points")
		}

		CheckBox
		{
			name: "showSigma"
			label: qsTr("State-noise standard deviation")
			checked: true
		}

		CheckBox
		{
			name: "showBeta"
			label: qsTr("Covariate coefficients")
			visible: hasCovariates
			info: qsTr("Include a named row for each covariate coefficient in the overview table.")
		}
	}
}
