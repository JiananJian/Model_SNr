#include "io_data.h"

int reset() {
	return std::remove(neuron_file) | std::remove(output_file);
}

void os_data(std::ostream& os, int N, int n, double** data) {
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < n; i++)
            os << data[i][j] << '\t';
        os << std::endl;
    }
}

void print_data(int N, int n, double** data) {
    os_data(std::cout, N, n, data);
}

void os_data(std::ostream& os, std::vector<double> data) {
    for (int j = 0; j < data.size(); j++) {
        os << data[j] << '\t';
        os << std::endl;
    }
}

int write_data(int N, int n, double** data) {
	std::ofstream outputFile(output_file, std::ios::app);
    if (outputFile.is_open()) {
        os_data(outputFile, N, n, data);
        outputFile.close();
        return 0;
    } else return -1;
}

void show_data(int N, double* x, double* y) {
    Gnuplot gp;
    std::vector<double> X(x, x + N);
    std::vector<double> Y(y, y + N);
    gp << "plot '-' with lines" << std::endl;
    gp.send1d(boost::make_tuple(X, Y));
}

void show_data(Gnuplot& gp, const char* filename, int x, int y) {
    gp << "plot \"" << filename << "\" using " << x << ":" << y << " with lines" << std::endl;
}
void show_data(const char* filename, int x, int y) {
    Gnuplot gp;
    gp << "plot \"" << filename << "\" using " << x << ":" << y << " with lines" << std::endl;
}

int save(const Neuron& neuron, const char* filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        file.write(reinterpret_cast<const char*>(&neuron), sizeof(neuron));
        return 0;
    } else return -1;
}
int load(Neuron& neuron, const char* filename) {
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open()) {
        Stim* GPe = neuron.GPe;
        Stim* Str = neuron.Str;
        file.read(reinterpret_cast<char*>(&neuron), sizeof(neuron));
        neuron.GPe = GPe;
        neuron.Str = Str;
        return 0;
    } else return -1;
}


std::ostream& operator<<(std::ostream& os, const Neuron& obj) {
#define PRINTVAR(x) os << #x << " = " << obj.x << std::endl
    PRINTVAR(time);
    //PRINTVAR(g_KCC2);
    //PRINTVAR(g_tonic);
    PRINTVAR(Cl_som);
    PRINTVAR(E_GABA_som);
    PRINTVAR(g_GABA_som);
    PRINTVAR(Cl_som);
    PRINTVAR(V_s);
    double E_Cl_som = V_T * log(Cl_out / obj.Cl_som) / z_Cl;
    double E_HCO3 = V_T * log(HCO3_out / HCO3_in) / z_HCO3; // -20mV
    double chi_som = (E_HCO3 - obj.E_GABA_som) / (E_HCO3 - E_Cl_som);
    double d_Cl = 1e3 * obj.alpha_Cl_som * (-chi_som * (obj.g_GABA_som) * (obj.V_s - E_Cl_som));

    os << "E_Cl_som" << " = " << E_Cl_som << std::endl;
    os << "d_Cl" << " = " << d_Cl << std::endl;
    //PRINTVAR(E_GABA_den);
    //PRINTVAR(V_d);
    //PRINTVAR(m_Na_f);
    //PRINTVAR(h_Na_f);
    //PRINTVAR(s_Na_f);
    //PRINTVAR(m_Na_p);
    //PRINTVAR(h_Na_p);
    //PRINTVAR(m_K);
    //PRINTVAR(h_K);
    //PRINTVAR(m_Ca);
    //PRINTVAR(h_Ca);
    return os;
}
