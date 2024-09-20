# Improved fluxonium readout through dynamic flux pulsing
by Taryn V. Stefanski, Figen Yilmaz, Eugene Y. Huang, Martijn F.S. Zwanenburg, Siddharth Singh, Siyu Wang, Lukas J. Splitthoff, and Christian Kraglund Andersen

## Description
Relevant code used for data analysis and numerical modeling of flux-pulse-assisted readout. Brief explanations of files given below.

### .csv files:
- res_fit_results.csv : contains fitted quantities from readout resonator spectroscopy with qubit prepared in |0> or |1>. The data was fit using the resonator magnitude model provided by our proprietary library from Orange Quantum Systems. This .csv file contains the resonator frequencies, linewidths, and dispersive shifts as a function of DC current bias converted to external flux bias in units of magnetic flux quanta, and the resonator frequency when the qubit is prepared in |0> as a function of flux pulse amplitude. Relevant data folders: 20231110, 20231111, 20231112, 20231113.
- processed_data_fx8.csv : exported quantities produced by process_fluxonium_data.ipynb.
- RO.csv : Measured data presented in Fig. 3 showing assignment and SNR-limited fidelity as a function of integration time at sweet-spot and using flux-pulse-assisted readout with 400mV readout amplitude, equivalent to 52 resonator photons.

### .ipynb files:
- process_fluxonium_data.ipynb : imports fitted resonator data from res_fit_results.csv, extracts measured qubit frequency and relaxation time as a function of external flux bias, and uses scQubits to model a fluxonium coupled to a resonator in order to fit the data and extract Ej, Ec, El, g, and (bare) resonator frequency. These values, along with qubit frequency and T1, are exported to processed_data_fx8.csv. Relevant data folders: 20231110, 20231111, 20231112, 20231113.
- exp_dispersive_shift_fitted.ipynb : imports res_fit_results.csv and processed_data_fx8.csv. Extracted energy parameters used to model dispersive shift and qubit-resonator detunings using code from https://github.com/AndersenQubitLab/FluxPulseAssistedFluxonium. Resonator frequency, qubit frequency, and dispersive shift data plotted in comparison to model prediction.
- process_FPARO_data.ipynb : imports and plots readout fidelity versus flux pulse amplitude for 5 different readout powers (100, 150, 200, 300, 400 mV), readout fidelity versus readout frequency at sweet-spot, and readout fidelity as a function of integration time for 5 different powers with the same, and variable, flux pulse amplitudes. Data pertaining to readout as a function of integration time using the highest readout amplitude (400 mV) which is shown in Fig. 3 is exported to RO.csv. Relevant data folders: 20231115, 20231116, 20231117, 20231118, 20231119, 20231201. 
- readout_dynamics_numerics.ipynb : imports data from processed_data_fx8.csv, res_fit_results.csv, and RO.csv to compare measured data with theory predictions from numerically solving the Langevin equation to obtain the intracavity field and resulting output field, adapted from the work in 'Flux-pulse-assisted readout of a fluxonium qubit.' This generates the theory curves of Fig. 3.

- meas_efficiency.ipynb : imports, plots, and analyzes coherence data, i.e. variable strength single-shot readout and Ramsey versus measurement pulse amplitude. Subsequently extracts measurement efficiency and imports data from res_fit_results.csv to convert readout amplitude to average resonator photon number. Additionally imports res_fit_results.csv to provide conversion of flux pulse amplitude to effective external flux bias in units of magnetic flux quanta. Relevant data folders: 20231129

### .py files
- ROfunctions.py : contains functions from https://github.com/AndersenQubitLab/FluxPulseAssistedFluxonium for diagonalizing full Hamiltonian, truncating, and calculating dispersive shift and detunings. Used in exp_dispersive_shift_fitted.ipynb.
- SSRO_functions.py : contains functions for processing single shot readout data. Used in process_FPARO_data.ipynb and meas_efficiency.ipynb.

## Requirements
- [Python == 3.8.8]
- [Cython == 0.29.20]
- [NumPy == 1.23.4]
- [SciPy == 1.8.1]
- [QuTiP == 4.7.0]
- [scQubits == 3.1.0]
- [Quantify-core == 0.7.4]
- [Matplotlib == 3.7.3]
- [Pandas == 1.2.4]
