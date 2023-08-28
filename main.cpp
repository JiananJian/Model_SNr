#include "utilities.h"
#include <random>

int main() {
	Neuron neuron; 

	// basic_characteristics(neuron); // Fig 1
	// g_SK_vs_rate(neuron);
	// voltage_clamp(neuron, -50);
	// current_clamp(neuron);
	test_inhibition(neuron, "GPe"); // Fig 4,5
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

	return 0;
}

