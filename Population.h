#pragma once

#include "Neuron.h"

class Population {
public:
    int size; // number of neurons
    Neuron* neurons; // array of neurons
    bool* spikes; // keep track of the spikes

    int n; // number of synapses 
    int* sources;
    int* targets;
    double* weights; // in nS/pF = 1/ms

	Stim GPe;
	Stim Str;

    Population(int size);
    ~Population();

    double time() const;
    void set_time(double time);
    void set_param(size_t param, double val);
    void set_param(size_t param, double* val);
    void step(double dt = 0.025, double* I_app = nullptr);
};