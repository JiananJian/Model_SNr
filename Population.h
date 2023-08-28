#pragma once

#include "Neuron.h"

class Population {
public:
    int size; // number of neurons
    Neuron* neurons;
    bool* spikes; // keep track of the spikes

    int n; // number of synapses 
    int* sources;
    int* targets;
    double* weights;

	Stim GPe;
	Stim Str;
    Population(int size = 2);
    ~Population();
    void step(double dt = 0.025);
};