#include "Population.h"

Population::Population(int size) : size(size) {
    neurons = new Neuron[size];
    spikes = new bool[size];
    for (int i = 0; i < size; i++) {
        neurons[i].GPe = &GPe;
        neurons[i].Str = &Str;
        spikes[i] = false;
    }
    n = 0;
    sources = new int[n];
    targets = new int[n];
    weights = new double[n];
}

Population::~Population() {
    delete[] neurons;
    delete[] spikes;
    delete[] sources;
    delete[] targets;
    delete[] weights;
}

void Population::step(double dt) {
    double* sum = new double[size]; 
    for (int i = 0; i < size; i++) 
        sum[i] = 0;
    for (int j = 0; j < n; j++) 
        if (spikes[sources[j]])
            sum[targets[j]] += weights[j];
    for (int i = 0; i < size; i++) 
        spikes[i] = neurons[i].step(dt, 0, sum[i]);
    delete[] sum;
}