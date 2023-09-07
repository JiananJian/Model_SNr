#include "utilities.h"
#include "Population.h"

#include "Sim.h"
int main() {
	Neuron neuron;
	test_inhib_w("GPe");

	return 0;

	reset();
	test_inhib_KCC2("GPe");
	return 0;
	double g_HCN_som = 20;
	double g_HCN_den = 20;
	while (g_HCN_den >= 0){
		g_HCN_som = 20;
		while (g_HCN_som >= 0) {
			Neuron neuron;
			neuron.W_Str = .4;
			neuron.W_GPe = .4;
			neuron.g_HCN_som = g_HCN_som;
			neuron.g_HCN_den = g_HCN_den;
			test_HCN(neuron, "");
			test_HCN(neuron, "GPe");
			test_HCN(neuron, "Str");
			std::cout << std::endl;
			g_HCN_som -= 5;
		}
		g_HCN_den -= 5;
	} 
	Plot gp;
	return 0;
	//neuron.Cl_som = 20;
	//neuron.Cl_den = 20;
	//neuron.W_Str = .4 * 160;
	//neuron.W_GPe = .2 * 6;
	// voltage_clamp(neuron, -70);

	// basic_characteristics(neuron); // Fig 1
	// g_SK_vs_rate(neuron);

	// current_clamp(neuron);
	// test_inhibition(neuron, "GPe"); // Fig 4,5
	// test_population();
	// g_HCN_vs_rate(neuron);

	// test_inhibition(neuron, "Str");

	return 0;
}

