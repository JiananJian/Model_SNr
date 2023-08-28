#include "Simulation.h"

int Simulation::run_sec() {
	double T = 1000; // duration in ms
	int N = (int)(T / dt);
	int n = 0;
	if (var) while (var[n]) n++;
	reset();
	int k = neuron.run(dt, N, n, var, val, spk, isi);
	if (var) write_data(N, n, val);

	std::cout << "spk(ms)" <<'\t' << "isi(ms)" << std::endl;
	print_data(k, spk, isi);
	std::cout << std::endl;
    return 0;
}

void Simulation::reset() {
	if (isi) delete[] isi;
	if (spk) delete[] spk;
	if (val) {
		int i = 0;
		while (var[i]) 
			delete[] val[i++];
		delete[] val;
	}
}

Simulation::~Simulation() {
	reset();
}
