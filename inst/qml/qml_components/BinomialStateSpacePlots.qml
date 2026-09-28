import QtQuick
import JASP.Controls

Section
{
	title: qsTr("Plots")
	columns: 1
	property bool hasCovariates: false

	CheckBox
	{
		name: "statePlot"
		label: qsTr("Probability trajectory")
		checked: true
		info: qsTr("Posterior mean probability over time with a 95% equal-tailed credible interval. The plot uses a percentage scale and conditions on all supplied observations.")
	}

	CheckBox
	{
		name: "posteriorDistPlot"
		label: qsTr("Posterior distribution of final probability")
		info: qsTr("Posterior probability distribution at the final time point in the supplied series, on a percentage scale. This is not a forecast beyond the series.")
	}

	CheckBox
	{
		name: "plotBeta"
		label: qsTr("Posterior distributions of covariate coefficients")
		visible: hasCovariates
		info: qsTr("A separate posterior distribution plot for each covariate coefficient. Each represents the change in log-odds for a one-standard-deviation increase in that covariate, holding the other covariates and latent level fixed.")
	}
}
