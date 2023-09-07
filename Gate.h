#pragma once

#include <math.h>

/**
 * @brief Parameters and dynamics of a voltage-dependent gate.
 */
class Gate {
public:
	const double V_z; // in mV
	const double k_z; // in mV
	const double x_min; // dimensionless
	const double V_tau; // in mV
	const double tau_0; // in ms
	const double tau_1; // in ms
	const double sig_0; // in mV
	const double sig_1; // in mV

	double dz(double z, double V, double dt) const; // dimensionless
};