# MSc-Research-Project

This is a [RIVET](https://rivet.hepforge.org/) analysis written for the Dijet Photoproduction paper "High-Et Dijet Photoproduction at HERA" (inspireID: [753991](https://inspirehep.net/literature/753991)) to compare experimental data presented in this paper to simulations from monte carlo generator. This analysis includes the C++  analysis code, the information file, the plotting file, and the experimental data file (.yoda).

## The Analysis Code

The analysis is written using the RIVET namespace. It works by initializing histograms for each dataset, then checking kinematics of each collision event to ensure it meets the kinematic conditions mentioned in the paper which are as follows:
 - at least two jets
 - photon-proton center of mass energies in the range 142 < $W_{\gamma P}$ < 293
 - jet energies $E_T^{jet1}$ > 20 GeV and $E_T^{jet2}$ > 15 GeV
 - psuedorapadities of -1 < $\eta$ < 3 with at least one in the range -1 < $\eta$ < 2.5
The analysis then proceeds to calculate kinematic variables of interest including $\Delta \phi$, $\overline{E_T}$, $\overline{eta}$, $x_{\gamma}^{obs}$, and $x_p^{obs}$. Next, the analysis fills the each histograms with it's corresponding kinematic variable, and checks conditions if applied to that plot. Finally, the analysis scales each histogram. An example of plots output by this analysis being run with the [Pythia](https://pythia.org/) monte carlo generator is seen below.

![example-rivet-plots/ZEUS_2007_I753991/d01-x01-y01.png]
