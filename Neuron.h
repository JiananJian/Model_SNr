#pragma once

#include "Stim.h"
#include "Gate.h"

// Thermal voltage at 308K in mV
#define V_T 26.54 

// Ion charges in e
#define z_Ca 2.
#define z_Cl -1.  
#define z_HCO3 -1.
#define z_GABA -1.

// Permeabilities
#define p_Cl 4. // Cl- ion permeability
#define p_HCO3 1. // HCO3- ion permeability

// Ion channel channel delay constants in nS/pF = 1/(ms)
#define g_Na_f 35. // fast Na+ channel 
#define g_Na_p .175 // persistent Na+ channel 
#define g_K 50. // K+ channel 
#define g_Ca .7 // Ca++ channel 
#define g_leak .04  // leak channel
#define g_TRPC3 .1 // TRPC3 channel 

// Capacitance in pF
#define C_som 100.
#define C_den 40. 

// Reversal potentials in mV
#define E_Na 50. 
#define E_K -90. 
#define E_leak -60. 
#define E_TRPC3 -37. 
#define E_HCN -45. // or -30 maybe?

// Ion concentrations in mM
#define Cl_out 120. 
#define HCO3_in 11.8 
#define HCO3_out 25. 
#define Ca_out 4. 

// Other constants
#define k_SK .4e-3 // m_SK half-maximal activation in mM
#define n_SK 4. // m_SK Hill coefficient, dimensionless
#define Ca_min 5e-8 // in mM
#define tau_Ca 250. // time constant in ms
#define tau_SD 200. // time constant in ms
#define tau_DS 80. // time constant in ms
#define tau_GABA_som 3. // time constant in ms 
#define tau_GABA_den 7.2 // time constant in ms 
#define tau_D 1000. // time constant in ms
#define tau_F 1000. // time constant in ms
#define D_0 1. // dimensionless
#define F_0 .145 // dimensionless
#define alpha_D .565 // dimensionless
#define alpha_F .125 // dimensionless
#define D_min .67 // dimensionless
#define F_min 1. // dimensionless

class Neuron {
public:
	/* ******** state variables ******** */

	// simulation time in ms
	double time = 0;

	// Membrane potentials in mV
	double V_s = E_leak;
	double V_d = E_leak;

	// Gate activation levels, dimensionless
	double m_Na_f = 0.1;
	double h_Na_f = 0.9;
	double s_Na_f = 0.9;
	double m_Na_p = 0.01;
	double h_Na_p = 0.04;
	double m_K = 0.01;
	double h_K = 0.9;
	double m_Ca = 0.001;
	double h_Ca = 0.001;
	double m_HCN_som = 0.01;
	double m_HCN_den = 0.01;

	// Synapse plasticity, dimensionless
	double D = D_0;
	double F = F_0;

	// Ion concentrations in mM
	double Ca_in = 2.5e-4;
	double Cl_som = 6;
	double Cl_den = 6;

	// Ion channel channel delay constants in nS/pF = 1/ms
	double g_GABA_som = 0;
	double g_GABA_den = 0;

	/* ******** output variables ******** */
	double dVs_dt = 0; // in mV/ms
	double dVd_dt = 0; // in mV/ms
	double I_GABA_som = 0; // in mV/ms = pA/pF
	double I_GABA_den = 0; // in mV/ms = pA/pF
	double E_GABA_som = 0; // in mV
	double E_GABA_den = 0; // in mV

	/* ******** parameters ******** */

	// Stimulation parameters
	Stim* GPe = nullptr;
	Stim* Str = nullptr;

	// Ion channel channel delay constants in nS/pF = 1/(ms)
	double g_SK = 0.009; // Calcium-activated K+ channel, reduces the firing rate
	double g_KCC2 = 0; // between 0.0 to 0.4 nS/pF
	double g_tonic = 0; // between 0.0 to 1.0 nS/pF
	double g_HCN_som = 0;
	double g_HCN_den = 0;
	double W_Str = .4;
	double W_GPe = .2;
	double W_SNr = .1;

	// parameters
	double V_th = -35; // spike threshold in mV

	// Charge to concentration conversion factors in mM/fC
	double alpha_Ca = 0.925e-7; // 1e-8 in Phillips2020
	double alpha_Cl_som = 1.85e-7; // 1.77e-7 in Phillips2020
	double alpha_Cl_den = 2.3e-6; // 2.2125e-7 in Phillips2020 

	// Conductance in nS
	double g_C = 26.5; // dendrite-soma coupling

	// Methods    
	/**
     * @brief Neuron dynamics in a time step.
     * @param dt Time step size in ms.
	 * @param I_app Applied current in pA.
	 * @param SNr Inputs in nS/pF received from other SNr within this time step.
	 * @param GPe Whether GPe stim is on.
	 * @param Str Whether D1 stim is on.
     * @return Whether an AP generates.
     */
	bool step(double dt, double I_app, double SNr, bool GPe, bool Str);

	/**
     * @brief Neuron dynamics in a time step.
	 * The GPe and Str inputs are specified by the member objects this->GPe and this->Str
     * @param dt Time step size in ms.
	 * @param I_app Applied current in pA.
	 * @param SNr Inputs in nS/pF received from other SNr within this time step.
     * @return Whether an AP generates.
     */
	bool step(double dt = 0.025, double I_app = 0, double SNr = 0);

private:
	// Ion gate parameters
	const Gate prop_m_Na_f = {
		-30.2,	// V_z
		6.2,	// k_z
		0,		// x_min
		1,		// V_tau
		0.05,	// tau_0
		0.05,	// tau_1
		1,		// sig_0
		1,		// sig_1
	};
	const Gate prop_h_Na_f = {
		-63.3,	// V_z
		-8.1,	// k_z
		0,		// x_min
		-43.,	// V_tau
		0.59,	// tau_0
		35.1,	// tau_1
		10,		//sig_0
		-5,		//sig_1
	};
	const Gate prop_s_Na_f = {
		-30.,	// V_z
		-0.4,	// k_z
		0.15,		// x_min
		-40,	// V_tau
		10,	// tau_0
		50,	// tau_1
		18.3,		//sig_0
		-10,		//sig_1
	};	
	const Gate prop_m_Na_p = {
		-50,	// V_z
		3,	// k_z
		0,		// x_min
		-42.6,	// V_tau
		0.03,	// tau_0
		0.146,	// tau_1
		14.4,		//sig_0
		-14.4,		//sig_1
	};
	const Gate prop_h_Na_p = {
		-57,	// V_z
		-4,	// k_z
		0.154,		// x_min
		-34,	// V_tau
		10,	// tau_0
		17,	// tau_1
		26,		//sig_0
		-31.9,		//sig_1
	};
	const Gate prop_m_K = {
		-26,	// V_z
		7.8,	// k_z
		0,		// x_min
		-26,	// V_tau
		0.1,	// tau_0
		14,	// tau_1
		13,		//sig_0
		-12,		//sig_1
	};
	const Gate prop_h_K = {
		-20,	// V_z
		-10,	// k_z
		0.6,		// x_min
		0,	// V_tau
		5,	// tau_0
		20,	// tau_1
		10,		//sig_0
		-10,		//sig_1
	};
	const Gate prop_m_Ca = {
		-27.5,	// V_z
		3,	// k_z
		0,		// x_min
		0,	// V_tau
		0.5,	// tau_0
		0.5,	// tau_1
		1,		//sig_0
		1,		//sig_1
	};
	const Gate prop_h_Ca = {
		-52.5,	// V_z
		-5.2,	// k_z
		0,		// x_min
		0,	// V_tau
		18,	// tau_0
		18,	// tau_1
		1,		//sig_0
		1,		//sig_1
	};
	const Gate prop_m_HCN = {
		-76.4,	// V_z
		-3.3,	// k_z
		0,		// x_min
		-76.4,	// V_tau
		0,		// tau_0
		3625,	// tau_1
		6.56,		// sig_0
		-7.48,		// sig_1
	};
};


