import QtQuick
import JASP
import JASP.Controls

Section
{
	title: qsTr("Advanced")
	columns: 1

	Group
	{
		title: qsTr("MCMC")

		IntegerField
		{
			id: warmup
			name: "advancedMcmcBurnin"
			label: qsTr("Warmup iterations per chain")
			defaultValue: 2000
			min: 0
			max: 2147483645
			info: qsTr("Iterations used for sampler adaptation and discarded before posterior summaries. Sampling iterations are additional to warmup.")
		}

		IntegerField
		{
			id: samples
			name: "advancedMcmcSamples"
			label: qsTr("Sampling iterations per chain")
			defaultValue: 5000
			min: 2
			max: 2147483647 - warmup.value
			info: qsTr("Post-warmup iterations per chain, before thinning. At least two draws per chain must remain after thinning.")
		}

		IntegerField
		{
			name: "advancedMcmcChains"
			label: qsTr("Chains")
			defaultValue: 3
			min: 1
			info: qsTr("Number of independent Markov chains. Multiple chains are needed to assess between-chain convergence.")
		}

		IntegerField
		{
			name: "advancedMcmcThin"
			label: qsTr("Thinning")
			defaultValue: 1
			min: 1
			max: Math.floor(samples.value / 2)
			info: qsTr("Retain every specified iteration after warmup. One keeps all sampling iterations.")
		}

		IntegerField
		{
			name: "advancedMcmcSeed"
			label: qsTr("Seed")
			defaultValue: 1
			min: 1
			max: 2147483647
			info: qsTr("Random seed for reproducible sampling within the same software environment.")
		}

		DoubleField
		{
			name: "advancedMcmcAdaptDelta"
			label: qsTr("Target acceptance probability")
			defaultValue: 0.97
			min: 0
			max: 1
			inclusive: JASP.None
			info: qsTr("Target acceptance probability for Hamiltonian Monte Carlo. Increasing it can reduce divergent transitions but may slow sampling.")
		}

		IntegerField
		{
			name: "advancedMcmcMaxTreeDepth"
			label: qsTr("Maximum tree depth")
			defaultValue: 15
			min: 1
			max: 20
			info: qsTr("Maximum tree depth for the sampler. Larger values permit more work per iteration and can substantially increase computation time.")
		}
	}
}
