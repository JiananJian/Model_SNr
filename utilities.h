#pragma once

#include "Population.h"
#include "Plot.h"
#include <random>

void basic_characteristics(Neuron& neuron);

void g_SK_vs_rate(Neuron& neuron);

void g_HCN_vs_rate(Neuron& neuron);

void voltage_clamp(Neuron& neuron, double V0);

void current_clamp(Neuron& neuron);

void test_inhibition(Neuron& neuron, const char* target);

void rate_inhibition(Neuron& neuron, const char* target);

void test_HCN(Neuron& neuron, const char* target);

void test_population();

void test_inhib_w(const char* target);

void test_inhib_gton(const char* target);

void test_inhib_KCC2(const char* target);



std::vector<double> run(Neuron& neuron, double start, double stop, const std::vector<size_t>& var, std::vector<std::vector<double>>& val);

// 
