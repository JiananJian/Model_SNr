#pragma once

#include "gnuplot-iostream.h"
#include <iostream>
#include <fstream>
#include <cstdio>

#define neuron_file "neuron.bin" // store the latest neuron state
#define output_file "output.txt" // store the temporal evolution of selected variables for further processing

// vector arithmetic
template<typename T>
std::vector<T> operator- (const std::vector<T>& a) {
    int m = a.size();
    std::vector<T> c;
    for (int i = 0; i < m; i++) c.push_back(-a[i]);
    return c;
}

template<typename T>
std::vector<T> operator+ (const std::vector<T>& a, const std::vector<T>& b) {
    int m = a.size(); int n = b.size();
    std::vector<T> c;
    if (m == n) for (int i = 0; i < m; i++) c.push_back(a[i] + b[i]);
    return c;
}

template<typename T>
std::vector<T> operator- (const std::vector<T>& a, const std::vector<T>& b) {
    int m = a.size(); int n = b.size();
    std::vector<T> c;
    if (m == n) for (int i = 0; i < m; i++) c.push_back(a[i] - b[i]);
    return c;
}

template<typename T>
std::vector<T> operator* (const T& a, const std::vector<T>& b) {
    int n = b.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a * b[i]);
    return c;
}

template<typename T>
std::vector<T> operator+ (const std::vector<T>& a, const T& b) {
    int m = a.size();
    std::vector<T> c;
    for (int i = 0; i < m; i++) c.push_back(a[i] + b);
    return c;
}

template<typename T>
std::vector<T> operator- (const std::vector<T>& a, const T& b) {
    return a + (-b);
}

template<typename T>
std::vector<T> operator* (const std::vector<T>& a, const T& b) {
    int n = a.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a[i] * b);
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

// methods gnuplot
extern void show_data(Gnuplot& gp, const char* filename = output_file, int x = 1, int y = 2);

extern void show_data(int N, double* x, double* y);

extern void show_data(const char* filename = output_file, int x = 1, int y = 2);