#include "Gate.h"

double Gate::dzdt(double z, double V) const {
	double z_0 = (1. - x_min) / (1. + exp((V_z - V) / k_z));
	double tau = tau_0 + (tau_1 - tau_0) / (exp((V_tau - V) / sig_0) + exp((V_tau - V) / sig_1));
	return (z_0 - z) / tau; // if tau too small, it implies z=z_0 for the sake of stability
}