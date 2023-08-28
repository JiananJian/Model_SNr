#include "Neuron.h"

/**
 * @brief Calculate single-neuron firing rate given applied current.
 * @param I_app Applied current in pA.
 * @return Firing rate in Hz.
 * @example 
 * 	std::cout << "I(pA)" << '\t' << "r(Hz)" << std::endl;
 *	for (int i = 0; i > -40; i--)
 *		std::cout << i << '\t' << I2r(i) << std::endl;
 */
double I2r(double I_app) {
	Neuron neuron; neuron.g_C = 4; neuron.g_SK = 0.0025;
	double dt = 0.025; // time step size in ms
	int size = 15; // number of spikes
	double* tsp = new double[size]; // spike timestamps
	int i = 0;
	while (i < size) {
		if (neuron.step(dt, I_app)) 
			tsp[i++] = neuron.time;
		if (neuron.time > 15e3) { // if the firing rate is roughly less than 1Hz
			delete[] tsp;
			return -1; 
		}
	}
	i = 5; 
	double rate = 1e3 / (tsp[size - 1] - tsp[i - 1]) * (size - i); // firing rate in Hz
	delete[] tsp;
	return rate;
}

/**
* One second simulation of the neuron model.
* Displays spike timestamps and inter-spike intervals as output.
* Appends #output_file.
*/
void run_sec(Neuron& neuron, const double** var = nullptr) {
	double dt = 0.025; // time step size in ms
	double T = 1000; // duration in ms
	int N = (int)(T / dt);
	int n = 0;
	if (var) while (var[n]) n++;
	double* spk = nullptr;
	double* isi = nullptr;
	double** val = nullptr;
	int k = neuron.run(dt, N, n, var, val, spk, isi);
	if (var) write_data(N, n, val);
	for (int i = 0; i < n; i++)
		delete[] val[i];
	delete[] val;

	//	std::cout << "spk(ms)" <<'\t' << "isi(ms)" << std::endl;
	//	print_data(k, spk, isi);
	//	std::cout << std::endl;
	delete[] spk;
	delete[] isi;
}

// study Fig 4A
void params_Cl() {
	Neuron neuron;
//	Gnuplot gp;
	const double* var[] = { &neuron.time, &neuron.Cl_som, nullptr };
	std::vector<double> g_KCC2, g_tonic, Cl_som, E_GABA;
	neuron.g_KCC2 = 0.4;
	neuron.g_tonic = 1.0;
	neuron.Cl_som = 5;
	double a;
	for (neuron.g_KCC2 = 0.0; neuron.g_KCC2 < 0.41 && neuron.g_tonic < 1.01; neuron.g_KCC2 += 0.01) {
//	for (neuron.g_tonic = 0.0; neuron.g_KCC2 < 0.41 && neuron.g_tonic < 1.01; neuron.g_tonic += 0.01) {
		a = 0;
		while (abs(a - neuron.Cl_som) > 1e-5) {
			a = neuron.Cl_som;
			run_sec(neuron);
//			std::cout << neuron << std::endl;
		}
		g_KCC2.push_back(neuron.g_KCC2);
		g_tonic.push_back(neuron.g_tonic);
		Cl_som.push_back(neuron.Cl_som);
		E_GABA.push_back(neuron.E_GABA_som);
		reset(); run_sec(neuron, var);
//		std::cout << neuron << std::endl;
//		std::string file = "neuron_" + std::to_string(neuron.g_KCC2) + "_" + std::to_string(neuron.g_tonic) + ".bin";
//		neuron.save(file.c_str());

//		gp << "plot ";
//		gp << "'" << output_file << "' using 1:2 with lines, ";
//		gp << std::endl;
	}
	print_data(Cl_som.size(), g_KCC2, g_tonic, Cl_som, E_GABA);
}


#include <random>
#include "Plot.h"
#include "utilities.h"
int main() {
	Neuron neuron; 

	// basic_characteristics(neuron); // Fig 1
	// g_SK_vs_rate(neuron);
	// voltage_clamp(neuron, -50);
	current_clamp(neuron);
	return 0;
	// Voltage clamp vs Current clamp
	if (1) {
		double V0 = -50;

		double dt = 0.025; // time step size in ms
		Stim stim_Str = { 1e3, 2e3, 25, dt };
		neuron.Str = &stim_Str;
		Stim stim_GPe = { 3e3, 4e3, 25, dt };
		neuron.GPe = &stim_GPe;

		const double* var[] = { &neuron.time, 
			&neuron.V_s, &neuron.V_d,
			& neuron.E_GABA_som,& neuron.E_GABA_den,
			& neuron.I_GABA_som,& neuron.I_GABA_den,
			& neuron.g_GABA_som, & neuron.g_GABA_den, 
			nullptr };
		int n = 0; if (var) while (var[n]) n++;

		neuron.time = -10e3;
		while (neuron.time < 0) {
			neuron.step();
		}
		std::vector<double>* val = new std::vector<double>[n];
		while (neuron.time < 5e3) {
			neuron.step();
			for (int j = 0; j < n; j++) val[j].push_back(*var[j]);
		}

		std::vector<double> I_DS = neuron.g_C * (val[2] - val[1]);

		Plot gp1i;
		gp1i.plot(val[0], I_DS, "from dendrite to soma");
		gp1i.xlabel("t (ms)").ylabel("I (pA)").xrange(0, 5e3);
		gp1i.flush();


		Plot gp1i3;
		gp1i3.plot(val[0], -val[7], "g_{som}").plot(val[0], -val[8], "g_{den}");
		gp1i3.xlabel("t (ms)").ylabel("g_{GABA}").xrange(0, 5e3);
		gp1i3.flush();


		Plot gp1v;
		gp1v.plot(val[0], val[1], "V_s").plot(val[0], val[2], "V_d");
		gp1v.xlabel("t (ms)").ylabel("V (mV)").xrange(0, 5e3);
		gp1v.flush();
		

		neuron.time = -10e3;
		while (neuron.time < 0) {
			neuron.step();
			neuron.V_s = V0; // voltage clamp
		}
		for (int j = 0; j < n; j++)	val[j].clear();
		while (neuron.time < 5e3) {
			neuron.step();
			for (int j = 0; j < n; j++) val[j].push_back(*var[j]);
			neuron.V_s = V0; // voltage clamp
		}
		Plot gp1i2;
		gp1i2.plot(val[0], -val[5] * C_som, "I_{som}").plot(val[0], -val[6] * C_den, "I_{den}");
		gp1i2.xlabel("t (ms)").ylabel("I_{GABA} (pA)").xrange(0, 5e3);
		gp1i2.flush();

		Plot gp3;
		gp3.plot(val[0], val[1], "V_s").plot(val[0], val[2], "V_d").
			plot(val[0], val[3], "E_GABA_som").plot(val[0], val[4], "E_GABA_den");
		gp3.xlabel("t (ms)").ylabel("V (mV)").xrange(0, 5e3);
		gp3.flush(); // display and save

		I_DS = neuron.g_C * (val[2] - val[1]);

		for (int i = 1; i < val[1].size(); i++) {
			val[1][i] -= V0;
			val[1][i] *= C_som / dt;
		}

		Plot gp2;
		gp2.plot(val[0], val[1], "").key(false);
		gp2.xlabel("t (ms)").ylabel("IPSC (pA)").xrange(0, 5e3);
		gp2.flush();//.save("V-Str-GPe.png"); // display and save


		Plot gp2i;
		gp2i.plot(val[0], I_DS, "").key(false);
		gp2i.xlabel("t (ms)").ylabel("I_{DS} (pA)").xrange(0, 5e3);
		gp2i.flush(); //.save("V-Str-GPe.png"); // display and save

		return 0;
	}

	// Fig 4,5
	if (0) {
		reset();
		Neuron neuron;
		neuron.g_KCC2 = 0.;
		neuron.g_tonic = 0.0;
		neuron.Cl_som = 5;

		double a = 0;
		while (abs(a - neuron.Cl_som) > 1e-5) {
			a = neuron.Cl_som;
			run_sec(neuron);
		}
		neuron.time = 0;

		const char* target = "GPe";
		Stim stim = { 1e3, 2e3, 40, 20 };
		if (~strcmp(target, "GPe"))
			neuron.GPe = &stim;
		if (~strcmp(target, "Str"))
			neuron.Str = &stim;

		const double* var[] = { &neuron.time, &neuron.V_s, &neuron.E_GABA_som, &neuron.E_GABA_den, &neuron.g_GABA_som, &neuron.Cl_som, &neuron.Cl_den, nullptr };
		int n = 0; if (var) while (var[n]) n++;
		std::vector<double>* val = new std::vector<double>[n];
		while (neuron.time < 3e3) {
			neuron.step();
			for (int j = 0; j < n; j++)
				val[j].push_back(*var[j]);
		}

		Plot gp1;
		gp1.plot(val[0], val[1], "V_s").plot(val[0], val[2], "E_GABA_som");
		gp1.xlabel("time (ms)").ylabel("V_s (mV)");
		gp1.flush();

		Plot gp2;
		gp2.plot(val[0], val[5], "Cl_som").plot(val[0], val[6], "Cl_den");
		gp2.xlabel("time (ms)").ylabel("Cl (mM)");
		gp2.flush();

		Plot gp3;
		gp3.plot(val[0], val[4], "g_GABA");
		gp3.xlabel("time (ms)").ylabel("g_GABA");
		gp3.flush();

		Plot gp4;
		gp4.plot(val[0], val[2], "E_GABA_som").plot(val[0], val[3], "E_GABA_den");
		gp4.xlabel("time (ms)").ylabel("E_GABA_som (mV)");
		gp4.flush();
	}

	return 0;

	if (1) {
		// Create a random engine
		std::random_device rd;
		std::mt19937 gen(rd());
		std::normal_distribution<double> dist;
		double I = dist(gen) * 10;

		Neuron neuron;
		const double* var[] = {&neuron.time, &neuron.V_s, nullptr};
			double dt = 0.025; // time step size in ms
			double T = 1000; // duration in ms
			int N = (int)(T / dt);
			int n = 0;
			if (var) while (var[n]) n++;

		double** val = new double* [n];
		for (int j = 0; j < n; j++)
			val[j] = new double[N];

		int k = 0;
		double* spk0 = new double[N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < n; j++)
				val[j][i] = *var[j];
			if (neuron.step(dt, I))
				spk0[k++] = neuron.time;
		}

		show_data(N, val[0], val[1]);
		return 0;
	}


	double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
	double y[5] = {1.0, 4.0, 9.0, 16.0, 25.0};

	Gnuplot gp;
	gp << "plot '-' with lines title 'x vs y'\n";
	gp.send1d(boost::make_tuple(x, y));

	return 0;
}

