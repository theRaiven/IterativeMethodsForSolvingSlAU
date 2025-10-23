#include "Header.h"
#include <cmath>
using namespace std;
void MultiplyMatrixByVector(const Matrix& A, const real* x, real* f, int n, int m, int offsets[], int sizes[])
{
    for (int i = 0; i < n; i++)
    {
        f[i] = 0.0;
        for (int d = 0; d < 9; ++d)
        {
            int offset = offsets[d];
            int j = i + offset;

            if (j < 0 || j >= n) continue;

            int i_min = max(0, -offset);
            int idx = i - i_min;

            if (idx < 0 || idx >= sizes[d]) continue;

            real aij{ 0.0 };
            switch (d)
            {
            case 0: aij = A.d4l[idx]; break;
            case 1: aij = A.d3l[idx]; break;
            case 2: aij = A.d2l[idx]; break;
            case 3: aij = A.d1l[idx]; break;
            case 4: aij = A.d0[idx];  break;
            case 5: aij = A.d1u[idx]; break;
            case 6: aij = A.d2u[idx]; break;
            case 7: aij = A.d3u[idx]; break;
            case 8: aij = A.d4u[idx]; break;
            }
            f[i] += aij * x[j];
        }
    }
}
real CalculateRelativeResidual(const Matrix& A, const real* x, const real* f, int n, int m, int offsets[], int sizes[])
{
    real* Ax = new real[n];
    MultiplyMatrixByVector(A, x, Ax, n, m, offsets, sizes);

    real normR_sq = 0.0;
    real normF_sq = 0.0;

    for (int i = 0; i < n; i++)
    {
        real ri = f[i] - Ax[i];
        normR_sq += ri * ri;
        normF_sq += f[i] * f[i];
    }

    delete[] Ax;

    if (normF_sq == 0) return 0.0;
    return sqrt(normR_sq / normF_sq);
}
void PrintResults(real omega, int iter, real normR, const real* x, const real* xTrue, int n)
{
    cout << fixed << setprecision(2) << "w = " << omega;
    cout << ", Iter = " << iter;
    cout << scientific << setprecision(15) << ", normR = " << normR << endl;

    cout << "i\t x_i" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i] << endl;
    }

    cout << "i\t x* _i - x_i" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << xTrue[i] - x[i] << endl;
    }
    cout << endl << endl;
}
void IterativeMethod(const Matrix& A, real* f, real* x, real* xTrue, int n, int m, real omega, int maxIter, bool isJacobi, int offsets[], int sizes[])
{
    
    if (n <= 0 || m < 0)
    {
        cerr << "Ошибка: некорректные размеры n или m!" << endl;
        return;
    }
    if (omega <= 0.0 || omega > 2.0)
    {
        cerr << "Ошибка: параметр релаксации ω должен быть в диапазоне (0, 2]!" << endl;
        return;
    }
    real* xNew = new real[n];
    real* r = new real[n];

    real normR{ 0.0 };
    int iter{ 0 };
    for (iter = 1; iter < maxIter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            real sum{ 0.0 };

            for (int d = 0; d < 9; ++d)
            {
                int offset = offsets[d];
                int j = i + offset;

                if (j < 0 || j >= n) continue;

                int i_min = max(0, -offset);
                int idx = i - i_min;

                if (idx < 0 || idx >= sizes[d]) continue;

                real aij{ 0.0 };
                switch (d)
                {
                case 0: aij = A.d4l[idx]; break;
                case 1: aij = A.d3l[idx]; break;
                case 2: aij = A.d2l[idx]; break;
                case 3: aij = A.d1l[idx]; break;
                case 5: aij = A.d1u[idx]; break;
                case 6: aij = A.d2u[idx]; break;
                case 7: aij = A.d3u[idx]; break;
                case 8: aij = A.d4u[idx]; break;
                }

                sum += aij * x[j];
            }

            real x_next = (f[i] - sum) / A.d0[i];

            if (isJacobi)
                xNew[i] = x[i] + omega * (x_next - x[i]);
            else
                x[i] = x[i] + omega * (x_next - x[i]);
        }

        if (isJacobi)
        {
            for (int i = 0; i < n; i++) x[i] = xNew[i];
        }

        normR = CalculateRelativeResidual(A, x, f, n, m, offsets, sizes);

        if (normR < relativeEPS<real>())
        {
            cout << "Достигнута требуемая точность решения." << endl;
            break;
        }

        if (iter >= maxIter - 1)
        {
            cout << "Достигнуто максимальное число итераций." << endl;
        }
    }

    PrintResults(omega, iter, normR, x, xTrue, n);

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
    }
    delete[] xNew;
    delete[] r;
}
void BlockRelaxation(const Matrix& A, real* f, real* x, real* xTrue, int n, int m, real omega, int maxIter, int blockSize, int offsets[], int sizes[])
{
    if (n <= 0 || blockSize <= 0 || blockSize > n)
    {
        cerr << "Ошибка: некорректный размер блока!" << endl;
        return;
    }

    real* r = new real[n];
    real normR{ 0.0 };

    int iter = 0;
    for (iter = 1; iter <= maxIter; iter++)
    {
        for (int i0 = 0; i0 < n; i0 += blockSize)
        {
            int i1 = min(i0 + blockSize, n);

            for (int i = i0; i < i1; i++)
            {
                real Ax{ 0.0 };
                for (int d = 0; d < 9; ++d)
                {
                    int offset = offsets[d];
                    int j = i + offset;

                    if (j < 0 || j >= n) continue;

                    int i_min = max(0, -offset);
                    int idx = i - i_min;

                    if (idx < 0 || idx >= sizes[d]) continue;

                    real aij{ 0.0 };
                    switch (d)
                    {
                    case 0: aij = A.d4l[idx]; break;
                    case 1: aij = A.d3l[idx]; break;
                    case 2: aij = A.d2l[idx]; break;
                    case 3: aij = A.d1l[idx]; break;
                    case 4: aij = A.d0[idx];  break;
                    case 5: aij = A.d1u[idx]; break;
                    case 6: aij = A.d2u[idx]; break;
                    case 7: aij = A.d3u[idx]; break;
                    case 8: aij = A.d4u[idx]; break;
                    }
                    Ax += aij * x[j];
                }
                r[i] = f[i] - Ax;
            }

            for (int i = i0; i < i1; i++)
            {
                real delta = (omega * r[i]) / A.d0[i];
                x[i] += delta;
            }
        }

        normR = CalculateRelativeResidual(A, x, f, n, m, offsets, sizes);

        if (normR < relativeEPS<real>())
        {
            cout << "Достигнута требуемая точность решения." << endl;
            break;
        }
    }
    real finalNormR = CalculateRelativeResidual(A, x, f, n, m, offsets, sizes);
    PrintResults(omega, iter, finalNormR, x, xTrue, n);

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
    }
    delete[] r;
}


void WorkingWithIterMethods()
{
	int n, m, maxIter;
    real* xTrue = nullptr;
    Matrix A;
    real* x0 = nullptr;
    real* y = nullptr;
    int sizes[9];
    
    try
    {
        InputDiagMatrix(A, xTrue, x0, y, n, m, maxIter, "Matrix.txt");
        CalculatingDimensionDiagonals(sizes, n, m);
        int offsets[]{ -4 - m, -3 - m, -2 - m, -1, 0, 1, 2 + m, 3 + m, 4 + m };
        cout << endl << n << ' ' << m << ' ' << maxIter << endl;
        PrintMatrixStruct(A, n, m);
        PrintMatrix(xTrue, n);
        PrintMatrix(x0, n);
        PrintMatrix(y, n);

        int methodChoice;
        cout << "\nВыберите метод:\n";
        cout << "1 — Метод Якоби\n";
        cout << "2 — Метод Гаусса–Зейделя\n";
        cout << "3 — Метод Блочной релаксации\n";
        cout << "Ваш выбор: ";
        cin >> methodChoice;
        cout << endl;

        switch (methodChoice)
        {
        case 1:
            cout << "\n=== Метод Якоби ===\n";
            for (float w = 0.01; w < 1.31; w += 0.1)
            {
                IterativeMethod(A, y, x0, xTrue, n, m, w, maxIter, true, offsets, sizes);
            }
            break;

        case 2:
            cout << "\n=== Метод Гаусса–Зейделя ===\n";
            for (float w = 0.01; w < 1.70; w += 0.1)
            {
                IterativeMethod(A, y, x0, xTrue, n, m, w, maxIter, false, offsets, sizes);
            }
            break;
        case 3:
            cout << "\n=== Метод блочной релаксации ===\n";
            int blockSize;
            cout << "Введите размер блока: ";
            cin >> blockSize;
            for (float w = 0.01; w < 2; w += 0.1)
            {
                BlockRelaxation(A, y, x0, xTrue, n, m, w, maxIter, blockSize, offsets, sizes);
            }
            break;
        default:
            throw runtime_error("Неверный выбор метода. Введите 1 или 2.");
        }
    }
    catch (exception& e)
    {
        cerr << "Ошибка: " << e.what() << endl;
    }
}