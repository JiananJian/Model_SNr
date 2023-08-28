#include "utilities.h"

void basic_characteristics(Neuron& neuron) {
	neuron.time = -1e3; while (neuron.time < 0) neuron.step(); // initialize

	const double* var[] = { &neuron.time, &neuron.V_s, &neuron.V_d, &neuron.dVs_dt, &neuron.dVd_dt, nullptr };
	int n = 0; if (var) while (var[n]) n++;
	std::vector<double>* val = new std::vector<double>[n];
	double spk = nan("");
	std::cout << "spk(ms)" << '\t' << "isi(ms)" << std::endl;
	while (neuron.time < 1e3) {
		for (int j = 0; j < n; j++) val[j].push_back(*var[j]);
		if (neuron.step()) {
			std::cout << neuron.time << '\t' << neuron.time - spk << std::endl;
			spk = neuron.time;
		}
	}
	std::cout << std::endl;

	std::string file;

	// Fig 1D
	Plot gp1d;
	file = "V_m phase portrait.png";
	gp1d.plot(val[1], val[3], "V_s").plot(val[2], val[4], "V_d");
	gp1d.xlabel("V_m (mV)").ylabel("dV_m/dt (mV/ms)");
	gp1d.flush().save(file.c_str()); // display and save

	// Fig 1B
	Plot gp1b;
	file = "V_m temporal evolution.png";
	gp1b.plot(val[0], val[1], "V_s").plot(val[0], val[2], "V_d");
	gp1b.xlabel("t (ms)").ylabel("V_m (mV)").xrange(0, 1e3);
	gp1b.flush().save(file.c_str()); // display and save

	// Fig 1C
	double dt = 0.025; // time step size in ms
	int size = 15; // number of spikes to record
	int k = 5; // discard the first k spikes
	double* tsp = new double[size]; // spike timestamps
	double I = 0, r = 1;
	std::vector<double> x, y;
	std::cout << "I(pA)" << '\t' << "r(Hz)" << std::endl;
	while (r) {
		neuron.time = 0;
		int i = 0;
		while (i < size) {
			if (neuron.step(dt, I)) tsp[i++] = neuron.time;
			if (i == size) {
				r = (size - k) / (tsp[size - 1] - tsp[k - 1]) * 1e3;
				break;
			}
			if (neuron.time > size * 1e3) { // if the firing rate is roughly less than 1Hz
				r = 0;
				break;
			}
		}
		std::cout << I << '\t' << r << std::endl;
		x.push_back(I--);
		y.push_back(r);
	}
	delete[] tsp;

	Plot gp1c;
	file = "f-I curve.png";
	gp1c.plot(x, y, "", "linespoints");
	gp1c.xlabel("I_{App} (pA)").ylabel("Frequency (Hz)").yrange(0, INFINITY).key(false);
	gp1c.flush().save(file.c_str()); // display and save
}

void g_SK_vs_rate(Neuron& neuron) {
	neuron.g_SK = 0;

	int size = 15; // number of spikes to record
	int k = 5; // discard the first k spikes
	double* tsp = new double[size]; // spike timestamps
	double r = 1;
	std::vector<double> x, y;
	std::cout << "g_SK(nS/pF)" << '\t' << "r(Hz)" << std::endl;
	while (r) {
		neuron.time = 0;
		int i = 0;
		while (i < size) {
			if (neuron.step()) tsp[i++] = neuron.time;
			if (i == size) {
				r = (size - k) / (tsp[size - 1] - tsp[k - 1]) * 1e3;
				std::cout << neuron.g_SK << '\t' << r << std::endl;
				x.push_back(neuron.g_SK); y.push_back(r);
				if (r > 5) neuron.g_SK += 1e-3;
				else if (r > 2) neuron.g_SK += 1e-2;
				else neuron.g_SK += 1e-1;
				break;
			}
			if (neuron.time > size * 1e3) { // if the firing rate is roughly less than 1Hz
				r = 0;
				break;
			}
		}
	}
	delete[] tsp;

	Plot gp1;
	gp1.plot(x, y, "", "lines");
	gp1.xlabel("g_{SK} (nS/pF)").ylabel("Frequency (Hz)").yrange(0, INFINITY).key(false);
	gp1.flush(); // display 
}

void voltage_clamp(Neuron& neuron, double V0) {

	double dt = 0.025; 
	Stim stim_Str = { 1e3, 2e3, 20, 0.025 };
	neuron.Str = &stim_Str;
	Stim stim_GPe = { 3e3, 4e3, 20, 0.025 };
	neuron.GPe = &stim_GPe;

	const double* var[] = { &neuron.time,
		&neuron.V_s, &neuron.V_d,
		&neuron.E_GABA_som,&neuron.E_GABA_den,
		&neuron.I_GABA_som,&neuron.I_GABA_den,
		&neuron.g_GABA_som, &neuron.g_GABA_den,
		nullptr };
	int n = 0; if (var) while (var[n]) n++;
	std::vector<double>* val = new std::vector<double>[n];

	neuron.time = -10e3;
	while (neuron.time < 0) {
		neuron.step(dt);
		neuron.V_s = V0; // voltage clamp
	}
	while (neuron.time < 5e3) {
		neuron.step(dt);
		for (int j = 0; j < n; j++) val[j].push_back(*var[j]);
		neuron.V_s = V0; // voltage clamp
	}

	std::vector<double> I_DS = neuron.g_C * (val[2] - V0);
	std::vector<double> IPSC = (val[1] - V0) * (C_som / dt);
	std::vector<double> I_som = -val[5] * C_som; // I_GABA in pA
	std::vector<double> I_den = -val[6] * C_den; // I_GABA in pA

	Plot gp0;
	gp0.plot(val[0], V0, "V_s").plot(val[0], val[2], "V_d").
		plot(val[0], val[3], "E_GABA_som").plot(val[0], val[4], "E_GABA_den");
	gp0.xlabel("t (ms)").ylabel("V (mV)").xrange(0, 5e3);
	gp0.flush();

	Plot gp1;
	gp1.plot(val[0], I_som, "I_{som}").plot(val[0], I_den, "I_{den}");
	gp1.xlabel("t (ms)").ylabel("I_{GABA} (pA)").xrange(0, 5e3);
	gp1.flush();

	Plot gp2;
	gp2.plot(val[0], IPSC, "").key(false);
	gp2.xlabel("t (ms)").ylabel("IPSC (pA)").xrange(0, 5e3);
	gp2.flush();//.save("IPSC.png"); // display and save

	Plot gp3;
	gp3.plot(val[0], I_DS, "").key(false);
	gp3.xlabel("t (ms)").ylabel("I_{DS} (pA)").xrange(0, 5e3);
	gp3.flush(); //.save("I_DS.png"); // display and save

}

void current_clamp(Neuron& neuron) {
	double dt = 0.025;
	Stim stim_Str = { 1e3, 2e3, 20, 0.025 };
	neuron.Str = &stim_Str;
	Stim stim_GPe = { 3e3, 4e3, 20, 0.025 };
	neuron.GPe = &stim_GPe;

	const double* var[] = { &neuron.time,
		&neuron.V_s, &neuron.V_d,
		&neuron.E_GABA_som,&neuron.E_GABA_den,
		&neuron.I_GABA_som,&neuron.I_GABA_den,
		&neuron.g_GABA_som, &neuron.g_GABA_den,
		nullptr };
	int n = 0; if (var) while (var[n]) n++;
	std::vector<double>* val = new std::vector<double>[n];

	neuron.time = -10e3;
	while (neuron.time < 0) {
		neuron.step(dt);
	}
	while (neuron.time < 5e3) {
		neuron.step(dt);
		for (int j = 0; j < n; j++) val[j].push_back(*var[j]);
	}

	std::vector<double> I_DS = neuron.g_C * (val[2] - val[1]);
	std::vector<double> I_som = -val[5] * C_som; // I_GABA in pA
	std::vector<double> I_den = -val[6] * C_den; // I_GABA in pA

	Plot gp0;
	gp0.plot(val[0], val[1], "V_s").plot(val[0], val[2], "V_d").
		plot(val[0], val[3], "E_GABA_som").plot(val[0], val[4], "E_GABA_den");
	gp0.xlabel("t (ms)").ylabel("V (mV)").xrange(0, 5e3);
	gp0.flush();

	Plot gp1;
	gp1.plot(val[0], I_som, "I_{som}").plot(val[0], I_den, "I_{den}");
	gp1.xlabel("t (ms)").ylabel("I_{GABA} (pA)").xrange(0, 5e3);
	gp1.flush();


	Plot gp3;
	gp3.plot(val[0], I_DS, "").key(false);
	gp3.xlabel("t (ms)").ylabel("I_{DS} (pA)").xrange(0, 5e3);
	gp3.flush(); //.save("I_DS.png"); // display and save


}