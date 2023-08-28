#pragma once

#include "Neuron.h"
#include "Plot.h"

void basic_characteristics(Neuron& neuron);

void g_SK_vs_rate(Neuron& neuron);

void voltage_clamp(Neuron& neuron, double V0);

void current_clamp(Neuron& neuron);