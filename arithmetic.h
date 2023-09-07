#pragma once

#include<vector>

/* vector arithmetic */


template<typename T>
std::vector<T> operator- (const std::vector<T>& a) {
    int m = a.size();
    std::vector<T> c;
    for (int i = 0; i < m; i++) c.push_back(-a[i]);
    return c;
}

template<typename T>
std::vector<T> log(const std::vector<T>& a) {
    int m = a.size();
    std::vector<T> c;
    for (int i = 0; i < m; i++) c.push_back(log(a[i]));
    return c;
}

template<typename T>
void operator/= (std::vector<T>& a, const T& b) {
    int m = a.size();
    for (int i = 0; i < m; i++) a[i] /= b;
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
std::vector<T> operator+ (const T& a, const std::vector<T>& b) {
    int n = b.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a + b[i]);
    return c;
}

template<typename T>
std::vector<T> operator- (const T& a, const std::vector<T>& b) {
    int n = b.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a - b[i]);
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
std::vector<T> operator/ (const T& a, const std::vector<T>& b) {
    int n = b.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a / b[i]);
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
    int m = a.size();
    std::vector<T> c;
    for (int i = 0; i < m; i++) c.push_back(a[i] - b);
    return c;
}

template<typename T>
std::vector<T> operator* (const std::vector<T>& a, const T& b) {
    int n = a.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a[i] * b);
    return c;
}

template<typename T>
std::vector<T> operator/ (const std::vector<T>& a, const T& b) {
    int n = a.size();
    std::vector<T> c;
    for (int i = 0; i < n; i++) c.push_back(a[i] / b);
    return c;
}
