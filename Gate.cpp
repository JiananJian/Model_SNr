#include "Gate.h"

double Gate::dz(double z, double V, double dt) const {
	double z_0 = (1. - x_min) / (1. + exp((V_z - V) / k_z));
	double tau = tau_0 + (tau_1 - tau_0) / (exp((V_tau - V) / sig_0) + exp((V_tau - V) / sig_1));
	return (z_0 - z) * (1 - exp(-dt / tau)); // stable even when tau is small
	return (z_0 - z) / tau * dt; // unstable when tau is small 
}
