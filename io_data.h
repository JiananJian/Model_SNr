#pragma once

#include "gnuplot-iostream.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include "Neuron.h"

#define neuron_file "neuron.bin" // store the latest neuron state
#define output_file "output.txt" // store the temporal evolution of selected variables for further processing

template<typename T>
std::vector<T> sampling (const std::vector<T>& a, int n = 10000) {
    int m = a.size();
    if (m <= n) return a;
    std::vector<T> c;
    for (int i = 0; i < m; i += m / n) c.push_back(a[i]);
    return c;
}


/**
 * @brief Delete #neuron_file and #output_file.
 */
extern int reset();

/**
 * @brief Stream multiple 1D arrays.
 * @param os Output stream object. 
 * @param N Number of rows of outputs.
 * @param args Arrays of data. 
 */
template<typename... Args>
void os_data(std::ostream& os, int N, Args... args) {
	for (int i = 0; i < N; i++) {
		((os << args[i] << '\t'), ...);
		os << std::endl;
	}
}

template<typename... Args>
int write_data(int N, Args... args) {
	std::ofstream outputFile(output_file, std::ios::app);
    if (outputFile.is_open()) {
        os_data(outputFile, N, args...);
        outputFile.close();
        return 0;
    } else {
        return -1;
    }
}

template<typename... Args>
void print_data(int N, Args... args) {
    os_data(std::cout, N, args...);
}

/**
 * @brief Stream a 2D arrays.
 * @param os Output stream object. 
 * @param N Number of rows of outputs.
 * @param n Number of columns of outputs.
 * @param data A 2D array of data. 
 */
extern void os_data(std::ostream& os, int N, int n, double** data);

extern void print_data(int N, int n, double** data);

extern int write_data(int N, int n, double** data);


// vectors
template<typename T, typename... Args>
void os_data(std::ostream& os, const std::vector<T>& arg0, const Args&... args) {
    for (int i = 0; i < arg0.size(); i++) {
        os << arg0[i] << '\t';
        ((os << args[i] << '\t'), ...);
        os << std::endl;
    }
}

template<typename T, typename... Args>
void print_data(const std::vector<T>& arg0, const Args&... args) {
    os_data(std::cout, arg0, args...);
}

template<typename T, typename... Args>
int write_data(const std::vector<T>& arg0, const Args&... args) {
    std::ofstream outputFile(output_file, std::ios::app);
    if (outputFile.is_open()) {
        os_data(outputFile, arg0, args...);
        outputFile.close();
        return 0;
    } else return -1;
}

template<typename T>
void os_data(std::ostream& os, const std::vector<T>* const args, const int n) {
    for (int i = 0; i < args->size(); i++) {
        for (int j = 0; j < n; j++)
            os << args[j][i] << '\t';
        os << std::endl;
    }
}

template<typename T>
void print_data(const std::vector<T>* const args, const int n) {
    os_data(std::cout, args, n);
}

template<typename T>
int write_data(const std::vector<T>* const args, const int n) {
    std::ofstream outputFile(output_file, std::ios::app);
    if (outputFile.is_open()) {
        os_data(outputFile, args, n);
        outputFile.close();
        return 0;
    } else return -1;
}

// methods gnuplot
extern void show_data(Gnuplot& gp, const char* filename = output_file, int x = 1, int y = 2);

extern void show_data(int N, double* x, double* y);

extern void show_data(const char* filename = output_file, int x = 1, int y = 2);

/**
 * @brief Save the neuron state into a binary file.
 * This method is sensitive to any modification to the member variables of this class.
 * The stimulation infomation cannot be saved or loaded.
 */
int save(const Neuron& neuron, const char* filename = neuron_file);

/**
 * @brief Load the neuron state from a binary file.
 * This method is sensitive to any modification to the member variables of this class.
 * The stimulation infomation cannot be saved or loaded.
 */
int load(Neuron& neuron, const char* filename = neuron_file);

/**
 * @brief Stream the values of selected state variables. Intended for debugging purposes.
 */
 std::ostream& operator<<(std::ostream& os, const Neuron& obj);
