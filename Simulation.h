#pragma once

#include "Neuron.h"
// do i want to save data? plot data? print data? for which pair of variables? or nothing? return a stream 
class Simulation {
public:
	Neuron neuron;
	const double** var = nullptr;
	double** val = nullptr;
	double* spk = nullptr;
	double* isi = nullptr;

	double dt = 0.025; // time step size in ms

	int run_sec();
	void reset();
	~Simulation();
};
