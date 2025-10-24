// Header.h

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

using real = double;
using realS = double;
const int accuracy = 15;

struct Matrix
{
    real* d4l = nullptr; // -4-m
    real* d3l = nullptr; // -3-m
    real* d2l = nullptr; // -2-m
    real* d1l = nullptr; // -1
    real* d0 = nullptr; // 0
    real* d1u = nullptr; // 1
    real* d2u = nullptr; // 2+m
    real* d3u = nullptr; // 3+m
    real* d4u = nullptr; // 4+m

    ~Matrix()
    {
        delete[] d4l;
        delete[] d3l;
        delete[] d2l;
        delete[] d1l;
        delete[] d0;
        delete[] d1u;
        delete[] d2u;
        delete[] d3u;
        delete[] d4u;
    }
};

//
template<typename T>
constexpr T relativeEPS();

template<>
constexpr float relativeEPS<float>() { return 1e-15f; }

template<>
constexpr double relativeEPS<double>() { return 1e-15; }


void InputDiagMatrix(Matrix& A, real*& xTrue, real*& x0, real*& y, int& n, int& m, string filename);
void CalculatingDimensionDiagonals(int sizes[], int n, int m);
void PrintMatrixStruct(const Matrix& A, int n, int m);

void PrintMatrix(real* A, int n);

void WorkingWithIterMethods();
