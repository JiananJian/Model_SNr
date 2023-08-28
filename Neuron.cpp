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
	double I_HCN = g_HCN * m_HCN * (V_s - E_HCN);

	double m_SK = 1. / (1. + pow(k_SK / Ca_in, n_SK));
	double I_SK = g_SK * m_SK * (V_s - E_K); // Calcium-activated K+ current 

	double I_TRPC3 = g_TRPC3 * (V_d - E_TRPC3);
	double I_SD = g_C / C_den * (V_d - V_s);
	I_GABA_den = g_GABA_den * (V_d - E_GABA_den);

	// State variable updates
	time += dt;

	m_Na_f += dt * prop_m_Na_f.dzdt(m_Na_f, V_s);
	h_Na_f += dt * prop_h_Na_f.dzdt(h_Na_f, V_s);
	s_Na_f += dt * prop_s_Na_f.dzdt(s_Na_f, V_s);
	m_Na_p += dt * prop_m_Na_p.dzdt(m_Na_p, V_s);
	h_Na_p += dt * prop_h_Na_p.dzdt(h_Na_p, V_s);
	m_K += dt * prop_m_K.dzdt(m_K, V_s);
	h_K += dt * prop_h_K.dzdt(h_K, V_s);
	m_Ca += dt * prop_m_Ca.dzdt(m_Ca, V_s);
	h_Ca += dt * prop_h_Ca.dzdt(h_Ca, V_s);
	m_HCN += dt * prop_m_HCN.dzdt(m_HCN, V_s);

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
	dVs_dt = -(I_Na_f + I_Na_p + I_K + I_Ca + I_SK + I_leak + I_GABA_som + I_DS + I_HCN) + I_app / C_som;
	dVd_dt = -(I_TRPC3 + I_GABA_den + I_SD);
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

int Neuron::run(double dt, int N, int n, const double** var, double**& val, double*& spk, double*& isi) {
		val = new double* [n];
		for (int j = 0; j < n; j++)
			val[j] = new double[N];

		int k = 0;
		double* spk0 = new double[N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < n; j++)
				val[j][i] = *var[j];
			if (step(dt))
				spk0[k++] = time;
		}

		spk = new double[k];
		for (int i = 0; i < k; i++)
			spk[i] = spk0[i];
		delete[] spk0;

		isi = new double[k]; 
		if (k) isi[0] = 0;
		for (int i = 1; i < k; i++)
			isi[i] = spk[i] - spk[i - 1];

		return k;
}

int Neuron::save(const char *filename) const {
	std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
		file.write(reinterpret_cast<const char*>(this), sizeof(*this));
        return 0;
    } else return -1;
}

int Neuron::load(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open()) {
        file.read(reinterpret_cast<char*>(this), sizeof(*this));
		GPe = nullptr;
		Str = nullptr;
		return 0;
    } else return -1;
}

std::ostream& operator<<(std::ostream& os, const Neuron& obj) {
	#define PRINTVAR(x) os << #x << " = " << obj.x << std::endl
	PRINTVAR(time);
	//PRINTVAR(g_KCC2);
	//PRINTVAR(g_tonic);
	PRINTVAR(Cl_som);
	PRINTVAR(E_GABA_som);
	PRINTVAR(g_GABA_som);
	PRINTVAR(Cl_som);
	PRINTVAR(V_s);
	double E_Cl_som = V_T * log(Cl_out / obj.Cl_som) / z_Cl;
	double E_HCO3 = V_T * log(HCO3_out / HCO3_in) / z_HCO3; // -20mV
	double chi_som = (E_HCO3 - obj.E_GABA_som) / (E_HCO3 - E_Cl_som);
	double d_Cl = 1e3 * alpha_Cl_som * (-chi_som * (obj.g_GABA_som) * (obj.V_s - E_Cl_som));

	os << "E_Cl_som" << " = " << E_Cl_som << std::endl;
	os << "d_Cl" << " = " << d_Cl << std::endl;
	//PRINTVAR(E_GABA_den);
	//PRINTVAR(V_d);
	//PRINTVAR(m_Na_f);
	//PRINTVAR(h_Na_f);
	//PRINTVAR(s_Na_f);
	//PRINTVAR(m_Na_p);
	//PRINTVAR(h_Na_p);
	//PRINTVAR(m_K);
	//PRINTVAR(h_K);
	//PRINTVAR(m_Ca);
	//PRINTVAR(h_Ca);
    return os;
}
