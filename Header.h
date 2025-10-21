// Header.h

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

using real = double;
using realS = double;
const int accuracy = 15;

//
template<typename T>
constexpr T relativeEPS();

template<>
constexpr float relativeEPS<float>() { return 1e-15f; }

template<>
constexpr double relativeEPS<double>() { return 1e-15; }


void InputDiagMatrix(real**& A, real*& xTrue, real*& x0, real*& y, int& n, int& m, int& maxIter, string filename);
void CalculatingDimensionDiagonals(int sizes[], int n, int m);
template<class T>
void PrintMatrix(T** A, int n, int m)
{
    cout << n << endl;
    int sizes[9];

    int k{ 0 };
    for (int i = 3; i > 0; i--)
    {
        sizes[k++] = n - m - (i + 1);
    }


    sizes[k++] = n - 1;
    sizes[k++] = n;
    sizes[k++] = n - 1;

    for (int i = 1; i <= 3; i++)
    {
        sizes[k++] = n - m - (i + 1);
    }
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < sizes[i]; j++)
        {
            cout <<  A[i][j] << ' ';
        }
        cout << endl;
    }
}
template<class T>
void PrintMatrix(T* A, int n)
{
    for (int j = 0; j < n; j++)
    {
        cout << setprecision(accuracy) <<  A[j] << ' ';
    }
    cout << endl;

}

void IterativeMethod(real** A, real* f, real* x, real* xTrue, int n, int m, real omega, int maxIter, bool isJacobi);
void WorkingWithIterMethods();