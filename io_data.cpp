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
