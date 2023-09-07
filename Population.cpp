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

double Population::time() const {
    if (size)
        return neurons[0].time;
    else
        return 0;
}

void Population::set_time(double time) {
    for (int i = 0; i < size; i++)
        neurons[i].time = time;
}

void Population::set_param(size_t param, double val) {
    for (int i = 0; i < size; i++) {
        double* param_ptr = (double*)(&neurons[i]) + param / sizeof(double);
        *param_ptr = val;
    }
}

void Population::set_param(size_t param, double* val) {
    for (int i = 0; i < size; i++) {
        double* param_ptr = (double*)(&neurons[i]) + param / sizeof(double);
        *param_ptr = val[i];
    }
}

void Population::step(double dt, double* I_app) {
    double* sum = new double[size]; 
    for (int i = 0; i < size; i++) 
        sum[i] = 0;
    for (int j = 0; j < n; j++) 
        if (spikes[sources[j]])
            sum[targets[j]] += weights[j];
    for (int i = 0; i < size; i++) {
        if (I_app)
            spikes[i] = neurons[i].step(dt, I_app[i], sum[i]);
        else
            spikes[i] = neurons[i].step(dt, 0, sum[i]);
    }
    delete[] sum;
}
