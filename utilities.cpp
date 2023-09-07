#include "utilities.h"
#include "arithmetic.h"

void basic_characteristics(Neuron& neuron) {
	neuron.time = -1e3; while (neuron.time < 0) neuron.step(); // initialize

	const double* var[] = { &neuron.time, &neuron.V_s, &neuron.V_d, &neuron.dVs_dt, &neuron.dVd_dt, nullptr };
	int n = 0; if (var) while (var[n]) n++;	std::vector<double>* val = new std::vector<double>[n];
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
	gp1.flush().save("g_SK.png"); // display 
}

void g_HCN_vs_rate(Neuron& neuron) {
	neuron.g_HCN_som = 0;

	int size = 15; // number of spikes to record
	int k = 5; // discard the first k spikes
	double* tsp = new double[size]; // spike timestamps
	double r = 1;
	std::vector<double> x, y;
	std::cout << "g_HCN(nS/pF)" << '\t' << "r(Hz)" << std::endl;
	while (r < 50) {
		neuron.time = 0;
		int i = 0;
		while (i < size) {
			if (neuron.step()) tsp[i++] = neuron.time;
			if (i == size) {
				r = (size - k) / (tsp[size - 1] - tsp[k - 1]) * 1e3;
				std::cout << neuron.g_HCN_som << '\t' << r << std::endl;
				x.push_back(neuron.g_HCN_som); y.push_back(r);
				neuron.g_HCN_som += 1;
				break;
			}
		}
	}
	delete[] tsp;

	Plot gp1;
	gp1.plot(x, y, "", "lines");
	gp1.xlabel("g_{HCN} (nS/pF)").ylabel("Frequency (Hz)").yrange(0, INFINITY).key(false);
	gp1.flush().save("g_HCN.png"); // display 
}

void voltage_clamp(Neuron& neuron, double V0) {
	neuron.Cl_som = 20; // E_GABA = -45mV
	neuron.Cl_den = 20; // E_GABA = -45mV

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
	std::vector<double> IPSC = -(val[1] - V0) * (C_som / dt);
	std::vector<double> I_som = -val[5] * C_som; // I_GABA in pA
	std::vector<double> I_den = -val[6] * C_den; // I_GABA in pA

	Plot gp0;
	gp0.plot(val[0], V0, "V_s").plot(val[0], val[2], "V_d").
		plot(val[0], val[3], "E_{GABA}_{som}").plot(val[0], val[4], "E_{GABA}_{den}");
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
		plot(val[0], val[3], "E_{GABA}_{som}").plot(val[0], val[4], "E_{GABA}_{den}");
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

void test_inhibition(Neuron& neuron, const char* target) {

	neuron.g_KCC2 = 0.0;
	neuron.g_tonic = 0.0;
	neuron.Cl_som = 6;

	// stabilize Cl_som for given g_KCC2 and g_tonic
	double Cl_som = 0;
	while (abs(Cl_som - neuron.Cl_som) > 1e-5) {
		Cl_som = neuron.Cl_som;
		neuron.step();
	}

	neuron.time = 0;

	Stim stim = { 1e3, 2e3, 20, 0.025 };
	if (strcmp(target, "GPe") == 0) {
		neuron.GPe = &stim;
		neuron.Str = nullptr;
	}
	if (strcmp(target, "Str") == 0) {
		neuron.Str = &stim;
		neuron.GPe = nullptr;
	}

	const double* var[] = { &neuron.time, &neuron.V_s, 
		&neuron.E_GABA_som, &neuron.E_GABA_den, &neuron.g_GABA_som, 
		&neuron.Cl_som, &neuron.Cl_den, nullptr };
	int n = 0; if (var) while (var[n]) n++;	
	std::vector<double>* val = new std::vector<double>[n];

	while (neuron.time < 3e3) {
		neuron.step();
		for (int j = 0; j < n; j++)	val[j].push_back(*var[j]);
	}

	Plot gp1;
	gp1.plot(val[0], val[1], "V_s").plot(val[0], val[2], "E_{GABA}_{som}");
	gp1.xlabel("time (ms)").ylabel("V_s (mV)");
	gp1.flush();

	Plot gp3;
	gp3.plot(val[0], val[4], "g_{GABA}");
	gp3.xlabel("time (ms)").ylabel("g_{GABA}");
	gp3.flush();

	Plot gp2;
	gp2.plot(val[0], val[5], "Cl_{som}").plot(val[0], val[6], "Cl_{den}");
	gp2.xlabel("time (ms)").ylabel("Cl (mM)");
	gp2.flush();

	Plot gp4;
	gp4.plot(val[0], val[2], "E_{GABA}_{som}").plot(val[0], val[3], "E_{GABA}_{den}");
	gp4.xlabel("time (ms)").ylabel("E_{GABA}_{som} (mV)");
	gp4.flush();
}

void rate_inhibition(Neuron& neuron, const char* target) {

	Stim stim = { 2e3, 12e3, 20, 0.025 };
	double* W; std::string w_str;
	if (strcmp(target, "GPe") == 0) {
		neuron.GPe = &stim;
		W = &neuron.W_GPe;
		w_str = "W_{GPe}";
	}
	if (strcmp(target, "Str") == 0) { 
		neuron.Str = &stim;
		W = &neuron.W_Str;
		w_str = "W_{Str}";
	}
	std::vector<double> spk, isi, t;

	Plot gp; 
	gp.xlabel("time (ms)").ylabel("firing rate (Hz)");

	*W = 0.0001;
	while (*W < 10) {
		*W *= 10;
		spk.clear(); isi.clear(); t.clear();
		neuron.time = -1e3;	while (neuron.time < 0) neuron.step();

		while (neuron.time < 15e3) {
			if (neuron.step()) {
				spk.push_back(neuron.time);
				std::cout << neuron.time << std::endl;
			}
		}

		for (int i = 1; i < spk.size(); i++) {
			t.push_back(spk[i - 1]);
			t.push_back(spk[i]);
			isi.push_back(spk[i] - spk[i - 1]);
			isi.push_back(spk[i] - spk[i - 1]);
			std::cout << isi[i] << std::endl;
		}

		std::string str = w_str + " = " + std::to_string(*W);
		gp.plot(t, 1e3 / isi, str.c_str());
		
	}

	gp.flush();

}

void test_HCN(Neuron& neuron, const char* target) {

	double g_HCN_som = neuron.g_HCN_som;
	double g_HCN_den = neuron.g_HCN_den;
	double ratio = g_HCN_den / g_HCN_som;
	neuron.g_HCN_som = 0;
	neuron.g_HCN_den = 0;
	
	int size = 15; // number of spikes to record
	while (neuron.g_HCN_som < g_HCN_som) {
		neuron.time = 0;
		int i = 0;
		while (i < size) {
			if (neuron.step()) i++;
			if (i == size) {
				neuron.g_HCN_som += 1;
				neuron.g_HCN_den += ratio;
				break;
			}
		}
	}
	neuron.g_HCN_som = g_HCN_som;
	neuron.g_HCN_den = g_HCN_den;
	

	Stim stim = { 2e3, 12e3, 20, 0.025 };
	if (strcmp(target, "GPe") == 0) {
		neuron.GPe = &stim;
		neuron.Str = nullptr;
	}
	if (strcmp(target, "Str") == 0) {
		neuron.Str = &stim;
		neuron.GPe = nullptr;
	}
		
	neuron.time = -10e3; while (neuron.time < 0) neuron.step(); 
	
	std::vector<double> spk, isi;

		const double* var[] = { &neuron.time, &neuron.V_s, &neuron.V_d, &neuron.dVs_dt, &neuron.dVd_dt, nullptr };
		int n = 0; if (var) while (var[n]) n++;
		std::vector<double>* val = new std::vector<double>[n];

		while (neuron.time < 15e3) {
			if (neuron.step()) {
				spk.push_back(neuron.time);
			}
			for (int j = 0; j < n; j++)
				val[j].push_back(*var[j]);
		}
		int count = 0;
		for (int i = 0; i < spk.size(); i++) {
			if (spk[i] > 2e3 && spk[i] < 12e3) count++;
		}
		double rate_opto = count / 10.0;

		count = 0;
		for (int i = 0; i < spk.size(); i++) {
			if (spk[i] < 2e3) count++;
		}
		double rate_rest = count / 2.0;
		std::cout << neuron.g_HCN_som << ", \t" << neuron.g_HCN_den << ", \t" << rate_rest << ", \t" << rate_opto << std::endl;


		Plot gp;
		//gp.xlabel("time (ms)").ylabel("firing rate (Hz)");
		//gp.xrange(0, 15);

		//gp.plot(t / 1e3, 1e3 / isi, "");
		//gp.flush();

		//Plot gp2;
		////gp2.plot(val[1], val[3], "").plot(val[2], val[4], "");
		//gp2.plot(val[0], val[1], "").plot(val[0], val[2], "");
		//gp2.flush();
}

void test_population() {

	// Create a random engine
	std::random_device rd;
	std::mt19937 gen(rd());
	std::lognormal_distribution<double> lognormal_dist;
	std::normal_distribution<double> normal_dist;

	int size = 5;
	Population population(size);
	population.GPe = Stim(1e3, 2e3, 20, 0.025);
	for (int i = 0; i < size; i++) {
		population.neurons[i].g_SK = lognormal_dist(gen);
	}
	double* I_app = new double[size];
	std::vector<double>* spikes = new std::vector<double>[size];
	std::vector<double>* V_s = new std::vector<double>[size];
	std::vector<double> time;

	while (population.time() < 10e3) {
		for (int i = 0; i < size; i++) {
			I_app[i] = normal_dist(gen) * 6e3;
		}
		population.step(0.025, I_app);
		time.push_back(population.time());
		for (int i = 0; i < population.size; i++) {
			V_s[i].push_back(population.neurons[i].V_s);
			if (population.spikes[i]) {
				spikes[i].push_back(population.neurons[i].time);
			}
		}

	}
	for (int i = 0; i < population.size; i++) {
		for (int j = 0; j < spikes[i].size(); j++)
			std::cout << spikes[i][j] << ", ";
		std::cout << std::endl << std::endl;
	}
	Plot gp1;
	for (int i = 0; i < population.size; i++) {
		gp1.plot(spikes[i], i, "", "points pt '*' lc rgb 'black'");
	}
	gp1.yrange(-1, population.size).flush();

	Plot gp2;
	for (int i = 0; i < population.size; i++) {
		gp2.plot(time, V_s[i], ""); break;
	}
	gp2.flush();

}

void test_inhib_w(const char* target) { // rate inhib v2
	std::vector<double> param_val = { 0 }; // {0, 1e-3, 1e-2, 1e-1};
	while (param_val.back() < 5) {
		param_val.push_back(param_val.back() + 0.1);
	}
	int size = param_val.size(); 

	Neuron neuron;
	neuron.V_th = -25;
	size_t param;
	const char* label;
	Stim stim = { 10e3, 20e3, 20, 0.025 };
	if (strcmp(target, "GPe") == 0) {
		label = "W_{GPe}";
		param = offsetof(Neuron, W_GPe);
		neuron.GPe = &stim;
		neuron.Str = nullptr;
	} else if (strcmp(target, "Str") == 0) {
		label = "W_{Str}";
		param = offsetof(Neuron, W_Str);
		neuron.GPe = nullptr;
		neuron.Str = &stim;
	} else {
		label = "W_{GPe}";
		param = offsetof(Neuron, W_GPe);
		neuron.Str = nullptr;
		neuron.GPe = nullptr;
	}
	double* param_ptr = (double*)(&neuron) + param / sizeof(double);
	std::vector<size_t> var = { offsetof(Neuron, time) ,offsetof(Neuron, V_s) };
	std::vector<std::vector<double>> val;
	std::vector<double> spk;

	std::cout << "without input" << "\t" << "with input" << "\t" << std::endl;
	std::cout << "rate(Hz)" << "\t" << "rate(Hz)" << "\t" << "W(nS/pF)" << std::endl;
	double rate0, rate1;
	Plot gp1;
	for (int i = 0; i < param_val.size(); i++) {
		*param_ptr = param_val[i];
		val.clear();
		spk = run(neuron, -10e3, 30e3, var, val);
		rate0 = 0;
		rate1 = 0;
		for (int j = 0; j < spk.size(); j++) {
			if (spk[j] < stim.start) rate0++;
			else if (spk[j] < stim.stop) rate1++;
		}
		rate0 /= 10;
		rate1 /= 10;
		std::cout << rate0 << "\t\t" << rate1 << "\t\t" << *param_ptr << std::endl;
		gp1.plot(spk * 1e-3, param_val[i], "", "points pt '+' lc rgb 'black'");
	}
	
	gp1.xlabel("time(sec)").ylabel("W(nS/pF)");
	gp1.flush();

	/*Plot gp2;
	gp2.plot(param_val, rate1);
	gp2.xlabel("W(nS/pF)").ylabel("FR(Hz)");
	gp2.flush();*/

	//write_data(val[1]);

}

void test_inhib_KCC2(const char* target) {
	std::vector<double> param_val = { .4, .2, .1, 0.05,0.008,0.007,0.004,.003,.002,.001,0 };
	//std::vector<double> param_val = {1.0 / 1024};
	//while (param_val.back() < 5) {
	//	param_val.push_back(param_val.back() * 2);
	//}
	int size = param_val.size(); Population population(size);
	Stim stim = { 10e3, 20e3, 20, 0.025 };
	if (strcmp(target, "GPe") == 0) {
		population.GPe = stim;
		population.Str = Stim();
	} else if (strcmp(target, "Str") == 0) {
		population.GPe = Stim();
		population.Str = stim;
	} else {
		population.Str = Stim();
		population.GPe = Stim();
	}

	size_t param = offsetof(Neuron, g_KCC2);
	double** param_ptr = new double* [size];
	for (int i = 0; i < size; i++) {
		param_ptr[i] = (double*)(&population.neurons[i]) + param / sizeof(double);
		*param_ptr[i] = param_val[i];
		population.neurons[i].V_th = -25;
		population.neurons[i].g_tonic = 1;
	}

	double alpha_Cl_som = population.neurons[0].alpha_Cl_som;
	double alpha_Cl_den = population.neurons[0].alpha_Cl_den;
	population.set_param(offsetof(Neuron, alpha_Cl_som), alpha_Cl_som * 1000);
	population.set_param(offsetof(Neuron, alpha_Cl_den), alpha_Cl_den * 1000);
	population.set_time(-10e3); while (population.time() < 0) population.step();
	population.set_param(offsetof(Neuron, alpha_Cl_som), alpha_Cl_som);
	population.set_param(offsetof(Neuron, alpha_Cl_den), alpha_Cl_den);

	Neuron* neuron = &population.neurons[0];
	const double* var[] = { &neuron->time, &neuron->V_s, &neuron->V_d, &neuron->dVs_dt, &neuron->dVd_dt, nullptr };
	int n = 0; if (var) while (var[n]) n++;
	std::vector<double>* val = new std::vector<double>[size + 1];
	std::vector<double>* spk = new std::vector<double>[size];
	while (population.time() < 30e3) {
		population.step();
		for (int i = 0; i < size; i++) {
			val[i].push_back(population.neurons[i].V_s);
			if (population.spikes[i])
				spk[i].push_back(population.time());
		}
		val[size].push_back(population.time());
	}

	std::cout << "without input" << "\t" << "with input" << "\t" << std::endl;
	std::cout << "rate(Hz)" << "\t" << "rate(Hz)" << "\t" << "g(nS/pF)" << "\t" 
		<< "Cl_som(mM)" << "\t" << "Cl_den(mM)" << "\t" << "E_GABA_som(mV)" << "\t" << "E_GABA_den(mV)" << "\t"
		<< std::endl;
	double* rate0 = new double[size];
	double* rate1 = new double[size];
	for (int i = 0; i < size; i++) {
		rate0[i] = 0;
		rate1[i] = 0;
		for (int j = 0; j < spk[i].size(); j++) {
			if (spk[i][j] < 10e3) rate0[i]++;
			else if (spk[i][j] < 20e3) rate1[i]++;
		}
		rate0[i] /= 10;
		rate1[i] /= 10;
		std::cout << rate0[i] << "\t\t" << rate1[i] << "\t\t" << *param_ptr[i] << "\t\t" 
			<< population.neurons[i].Cl_som << "\t\t" << population.neurons[i].Cl_den << "\t\t"
			<< population.neurons[i].E_GABA_som << "\t\t" << population.neurons[i].E_GABA_den << "\t\t"
			<< std::endl;
	}

	Plot gp1;
	for (int i = 0; i < size; i++) {
		gp1.plot(spk[i] * 1e-3, log(param_val[i]) / log(2), "", "points pt '+' lc rgb 'black'");
	}
	gp1.xlabel("time(sec)").ylabel("g_{KCC2}(nS/pF): 2 to the power of ");
	gp1.flush();

	Plot gp2;
	gp2.plot(log(param_val) / log(2), rate1);
	gp2.xlabel("g_{KCC2}(nS/pF): 2 to the power of").ylabel("FR(Hz)");
	gp2.flush();

	write_data(val, size+1);

	delete[] val;
	delete[] spk;
	delete[] param_ptr;
	delete[] rate0;
	delete[] rate1;
}


void test_inhib_gton(const char* target) {

	std::vector<double> param_val = { 1.0 / 1024 };
	while (param_val.back() < 5) {
		param_val.push_back(param_val.back() * 2);
	}
	int size = param_val.size(); Population population(size);
	Stim stim = { 10e3, 20e3, 20, 0.025 };
	if (strcmp(target, "GPe") == 0) {
		population.GPe = stim;
		population.Str = Stim();
	} else if (strcmp(target, "Str") == 0) {
		population.GPe = Stim();
		population.Str = stim;
	} else {
		population.Str = Stim();
		population.GPe = Stim();
	}

	size_t param = offsetof(Neuron, g_tonic);
	double** param_ptr = new double* [size];
	for (int i = 0; i < size; i++) {
		param_ptr[i] = (double*)(&population.neurons[i]) + param / sizeof(double);
		*param_ptr[i] = param_val[i];
		population.neurons[i].V_th = -25;
		population.neurons[i].g_KCC2 = 1;
	}

	double alpha_Cl_som = population.neurons[0].alpha_Cl_som;
	double alpha_Cl_den = population.neurons[0].alpha_Cl_den;
	population.set_param(offsetof(Neuron, alpha_Cl_som), alpha_Cl_som * 1000);
	population.set_param(offsetof(Neuron, alpha_Cl_den), alpha_Cl_den * 1000);
	population.set_time(-10e3); while (population.time() < 0) population.step();
	population.set_param(offsetof(Neuron, alpha_Cl_som), alpha_Cl_som);
	population.set_param(offsetof(Neuron, alpha_Cl_den), alpha_Cl_den);

	Neuron* neuron = &population.neurons[0];
	const double* var[] = { &neuron->time, &neuron->V_s, &neuron->V_d, &neuron->dVs_dt, &neuron->dVd_dt, nullptr };
	int n = 0; if (var) while (var[n]) n++;
	std::vector<double>* val = new std::vector<double>[size + 1];
	std::vector<double>* spk = new std::vector<double>[size];
	while (population.time() < 30e3) {
		population.step();
		for (int i = 0; i < size; i++) {
			val[i].push_back(population.neurons[i].V_s);
			if (population.spikes[i])
				spk[i].push_back(population.time());
		}
		val[size].push_back(population.time());
	}

	std::cout << "without input" << "\t" << "with input" << "\t" << std::endl;
	std::cout << "rate(Hz)" << "\t" << "rate(Hz)" << "\t" << "g(nS/pF)" << "\t"
		<< "Cl_som(mM)" << "\t" << "Cl_den(mM)" << "\t" << "E_GABA_som(mV)" << "\t" << "E_GABA_den(mV)" << "\t"
		<< std::endl;
	double* rate0 = new double[size];
	double* rate1 = new double[size];
	for (int i = 0; i < size; i++) {
		rate0[i] = 0;
		rate1[i] = 0;
		for (int j = 0; j < spk[i].size(); j++) {
			if (spk[i][j] < 10e3) rate0[i]++;
			else if (spk[i][j] < 20e3) rate1[i]++;
		}
		rate0[i] /= 10;
		rate1[i] /= 10;
		std::cout << rate0[i] << "\t\t" << rate1[i] << "\t\t" << *param_ptr[i] << "\t\t"
			<< population.neurons[i].Cl_som << "\t\t" << population.neurons[i].Cl_den << "\t\t"
			<< population.neurons[i].E_GABA_som << "\t\t" << population.neurons[i].E_GABA_den << "\t\t"
			<< std::endl;
	}

	Plot gp1;
	for (int i = 0; i < size; i++) {
		gp1.plot(spk[i] * 1e-3, log(param_val[i]) / log(2), "", "points pt '+' lc rgb 'black'");
	}
	gp1.xlabel("time(sec)").ylabel("g_{tonic}(nS/pF): 2 to the power of ");
	gp1.flush();

	Plot gp2;
	gp2.plot(log(param_val) / log(2), rate1);
	gp2.xlabel("g_{tonic}(nS/pF): 2 to the power of").ylabel("FR(Hz)");
	gp2.flush();

	write_data(val, size + 1);

	delete[] val;
	delete[] spk;
	delete[] param_ptr;
	delete[] rate0;
	delete[] rate1;
}



std::vector<double> run(Neuron& neuron, double start, double stop, 
	const std::vector<size_t>& var, std::vector<std::vector<double>>& val) {
	neuron.time = start; while (neuron.time < 0) neuron.step();

	std::vector<double> spk;
	val.resize(var.size());
	std::vector<double*> ptr(var.size());
	for (int i = 0; i < var.size(); i++) {
		ptr[i] = (double*)(&neuron) + var[i] / sizeof(double);
	}
	while (neuron.time < stop) {
		if (neuron.step()) spk.push_back(neuron.time);
		for (int i = 0; i < var.size(); i++) {
			val[i].push_back(*ptr[i]);
		}
	}
	return spk;
}

