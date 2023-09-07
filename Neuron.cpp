#include "Neuron.h"

bool Neuron::step(double dt, double I_app, double SNr, bool GPe, bool Str) {
    // Reversal potentials in mV
	double E_Ca = V_T * log(Ca_out / Ca_in) / z_Ca; 
	double E_HCO3 = V_T * log(HCO3_out / HCO3_in) / z_HCO3; // -20mV
	double E_Cl_som = V_T * log(Cl_out / Cl_som) / z_Cl; 
	double E_Cl_den = V_T * log(Cl_out / Cl_den) / z_Cl; 
	E_GABA_som = V_T * log((p_Cl * Cl_out + p_HCO3 * HCO3_out) / (p_Cl * Cl_som + p_HCO3 * HCO3_in)) / z_GABA;
	E_GABA_den = V_T * log((p_Cl * Cl_out + p_HCO3 * HCO3_out) / (p_Cl * Cl_den + p_HCO3 * HCO3_in)) / z_GABA;

	// Outward currents in mV/ms
	double I_Na_f = g_Na_f * pow(m_Na_f, 3) * h_Na_f * s_Na_f * (V_s - E_Na); // fast Na+ current
	double I_Na_p = g_Na_p * pow(m_Na_p, 3) * h_Na_p * (V_s - E_Na); // persistent Na+ current
	double I_K = g_K * pow(m_K, 4) * h_K * (V_s - E_K); // K+ current
	double I_Ca = g_Ca * m_Ca * h_Ca * (V_s - E_Ca); // Ca++ current
	double I_leak = g_leak * (V_s - E_leak);
	I_GABA_som = g_GABA_som * (V_s - E_GABA_som);
	double I_DS = g_C / C_som * (V_s - V_d);
	double I_HCN_som = g_HCN_som * m_HCN_som * (V_s - E_HCN);
	double I_HCN_den = g_HCN_den * m_HCN_den * (V_d - E_HCN);

	double m_SK = 1. / (1. + pow(k_SK / Ca_in, n_SK));
	double I_SK = g_SK * m_SK * (V_s - E_K); // Calcium-activated K+ current 

	double I_TRPC3 = g_TRPC3 * (V_d - E_TRPC3);
	double I_SD = g_C / C_den * (V_d - V_s);
	I_GABA_den = g_GABA_den * (V_d - E_GABA_den);

	// State variable updates
	time += dt;

	m_Na_f += prop_m_Na_f.dz(m_Na_f, V_s, dt);
	h_Na_f += prop_h_Na_f.dz(h_Na_f, V_s, dt);
	s_Na_f += prop_s_Na_f.dz(s_Na_f, V_s, dt);
	m_Na_p += prop_m_Na_p.dz(m_Na_p, V_s, dt);
	h_Na_p += prop_h_Na_p.dz(h_Na_p, V_s, dt);
	m_K += prop_m_K.dz(m_K, V_s, dt);
	h_K += prop_h_K.dz(h_K, V_s, dt);
	m_Ca += prop_m_Ca.dz(m_Ca, V_s, dt);
	h_Ca += prop_h_Ca.dz(h_Ca, V_s, dt);
	m_HCN_som += prop_m_HCN.dz(m_HCN_som, V_s, dt);
	m_HCN_den += prop_m_HCN.dz(m_HCN_den, V_d, dt);

	Ca_in -= dt * (alpha_Ca * C_som * I_Ca + (Ca_in - Ca_min) / tau_Ca);

	double chi_som = (E_HCO3 - E_GABA_som) / (E_HCO3 - E_Cl_som);
	double chi_den = (E_HCO3 - E_GABA_den) / (E_HCO3 - E_Cl_den);
	Cl_som -= dt * alpha_Cl_som * C_som * (g_KCC2 * (E_Cl_som - E_K) - chi_som * (g_GABA_som + g_tonic) * (V_s - E_Cl_som));
	Cl_den -= dt * alpha_Cl_den * C_den * (g_KCC2 * (E_Cl_den - E_K) - chi_den * (g_GABA_den + g_tonic) * (V_d - E_Cl_den));
	Cl_som -= dt * (Cl_som - Cl_den) / tau_SD;
	Cl_den -= dt * (Cl_den - Cl_som) / tau_DS;

	g_GABA_som *= exp(-dt / tau_GABA_som);
	g_GABA_den *= exp(-dt / tau_GABA_den);

	g_GABA_som += SNr;
	
	if (GPe) {
		g_GABA_som += W_GPe * D;
		D -= alpha_D * (D - D_min);
	}
	if (Str) {
		g_GABA_den += W_Str * F;
		F -= alpha_F * (F - F_min);
	}
	D += dt * (D_0 - D) / tau_D;
	F += dt * (F_0 - F) / tau_F;

	double V_0 = V_s;
	dVs_dt = -(I_Na_f + I_Na_p + I_K + I_Ca + I_SK + I_leak + I_GABA_som + I_DS + I_HCN_som) + I_app / C_som;
	dVd_dt = -(I_TRPC3 + I_GABA_den + I_SD + I_HCN_den);
	V_s += dt * dVs_dt;
	V_d += dt * dVd_dt;
	return (V_0 < V_th) && (V_s >= V_th);
}

bool Neuron::step(double dt, double I_app, double SNr) {
	bool GPe_0 = false;
	bool Str_0 = false;
	if (GPe) GPe_0 = GPe->pulse(time);
	if (Str) Str_0 = Str->pulse(time);
	return step(dt, I_app, SNr, GPe_0, Str_0);
}

