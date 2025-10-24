#include "Header.h"
#include <cmath>

struct SolverParams
{
    real omega{0};
    int maxIter{10000};
    bool isJacobi{false};
    int blockSize{0};
};

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

    real normR_sq{ 0.0 };
    real normF_sq{ 0.0 };

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

void IterationStep(const Matrix& A, const real* f, real*& x, real* xNew, int n, int m, const SolverParams& params, int offsets[], int sizes[])
{
    for (int i = 0; i < n; i++)
    {
        real sum{ 0.0 };

        for (int d = 0; d < 9; ++d)
        {
            if (d == 4) continue;

            int offset{ offsets[d] };
            int j{ i + offset };
            if (j < 0 || j >= n) continue;

            int i_min{ max(0, -offset) };
            int idx{ i - i_min };
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
            default: continue;
            }

            sum += aij * x[j];
        }

        real x_next = (f[i] - sum) / A.d0[i];
        if (params.isJacobi)
            xNew[i] = x[i] + params.omega * (x_next - x[i]);
        else
            x[i] = x[i] + params.omega * (x_next - x[i]);
    }
}
void IterativeMethod(const Matrix& A, real* f, real*& x, real* xTrue, int n, int m, SolverParams params, int offsets[], int sizes[], real* xNew, real* r)
{
    if (n <= 0 || m < 0)
    {
        cerr << "Ошибка: некорректные размеры n или m!" << endl;
        return;
    }
    if (params.omega <= 0.0 || params.omega > 2.0)
    {
        cerr << "Ошибка: параметр релаксации ω должен быть в диапазоне (0, 2]!" << endl;
        return;
    }
    int iter{ 0 };
    real normR{ CalculateRelativeResidual(A, x, f, n, m, offsets, sizes) };
    if (normR < relativeEPS<real>())
    {
        cout << "Начальное приближение уже является решением." << endl;
        PrintResults(params.omega, iter, normR, x, xTrue, n);
        return;
    }

    for (iter = 1; iter < params.maxIter; iter++)
    {
        IterationStep(A, f, x, xNew, n, m, params, offsets, sizes);

        if (params.isJacobi)
        {
            for (int i = 0; i < n; i++) x[i] = xNew[i];
        }

        normR = CalculateRelativeResidual(A, x, f, n, m, offsets, sizes);

        if (normR < relativeEPS<real>())
        {
            cout << "Достигнута требуемая точность решения." << endl;
            break;
        }
    }
    if (iter >= params.maxIter - 1)
    {
        cout << "Достигнуто максимальное число итераций." << endl;
    }

    PrintResults(params.omega, iter, normR, x, xTrue, n);
    WriteResultsToFile(params.isJacobi ? "Метод Якоби" : "Метод Гаусса–Зейделя",
        params.omega, iter, normR, x, xTrue, n, "Result.txt");
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        xNew[i] = 0;
        r[i] = 0;
    }
}
void BlockRelaxation(const Matrix& A, real* f, real*& x, real* xTrue, int n, int m, SolverParams params, int offsets[], int sizes[], real* r)
{
    if (n <= 0 || params.blockSize <= 0 || params.blockSize > n)
    {
        cerr << "Ошибка: некорректный размер блока!" << endl;
        return;
    }

    int iter{ 0 };
    real normR{ CalculateRelativeResidual(A, x, f, n, m, offsets, sizes) };
    if (normR < relativeEPS<real>())
    {
        cout << "Начальное приближение уже является решением." << endl;
        PrintResults(params.omega, iter, normR, x, xTrue, n);
        return;
    }

    for (iter = 1; iter <= params.maxIter; iter++)
    {
        for (int i0 = 0; i0 < n; i0 += params.blockSize)
        {
            int i1{ min(i0 + params.blockSize, n) };

            for (int i = i0; i < i1; i++)
            {
                real Ax{ 0.0 };
                for (int d = 0; d < 9; ++d)
                {
                    int offset{ offsets[d] };
                    int j{ i + offset };

                    if (j < 0 || j >= n) continue;

                    int i_min{ max(0, -offset) };
                    int idx{ i - i_min };

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
                    Ax += aij * ((j >= i0 && j < i) ? x[j] : x[j]);
                }
                r[i] = f[i] - Ax;

                if (A.d0[i] == 0) throw runtime_error("Деление на Нуль!");
                real delta = (params.omega * r[i]) / A.d0[i];
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
    PrintResults(params.omega, iter, finalNormR, x, xTrue, n);
    WriteResultsToFile("Метод блочной релаксации", params.omega, iter, finalNormR, x, xTrue, n, "Result.txt");

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        r[i] = 0;
    }
}


void WorkingWithIterMethods()
{
	int n, m;
    real* xTrue = nullptr;
    Matrix A;
    real* x0 = nullptr;
    real* y = nullptr;
    int sizes[9];
    
    try
    {
        InputDiagMatrix(A, xTrue, x0, y, n, m, "Matrix.txt");
        CalculatingDimensionDiagonals(sizes, n, m);
        SolverParams params;
        int offsets[]{ -4 - m, -3 - m, -2 - m, -1, 0, 1, 2 + m, 3 + m, 4 + m };
        real* xNew = new real[n];
        real* r = new real[n];
        cout << endl << n << ' ' << m << endl;
        PrintMatrixStruct(A, n, m);
        PrintMatrix(xTrue, n);
        PrintMatrix(x0, n);
        PrintMatrix(y, n);

        bool exit{true};
        while (exit)
        {
            int methodChoice;
            cout << "\nВыберите метод:\n";
            cout << "1 — Метод Якоби\n";
            cout << "2 — Метод Гаусса–Зейделя\n";
            cout << "3 — Метод Блочной релаксации\n";
            cout << "0 — Выход\n";
            cout << "Ваш выбор: ";
            cin >> methodChoice;
            cout << endl;


            float w0, w1, v;
        
            switch (methodChoice)
            {
            case 1:
                cout << "\n=== Метод Якоби ===\n";
                params.maxIter = 100000;
                params.isJacobi = true;
                cout << "Задайте диапазон omega через пробел(пример: 0.01 1.61): ";
                float w0, w1, v;
                cin >> w0 >> w1;
                cout << "Задайте шаг omega (пример: 0.01): ";
                cin >> v;
                for (params.omega = w0; params.omega < w1; params.omega += v)
                {
                    IterativeMethod(A, y, x0, xTrue, n, m, params, offsets, sizes, xNew, r);
                }
                break;

            case 2:
                cout << "\n=== Метод Гаусса–Зейделя ===\n";
                params.maxIter = 1000000;
                params.isJacobi = false;
                cout << "Задайте диапазон omega через пробел(пример: 0.01 1.61): ";
                cin >> w0 >> w1;
                cout << "Задайте шаг omega (пример: 0.01): ";
                cin >> v;
                for (params.omega = w0; params.omega < w1; params.omega += v)
                {
                    IterativeMethod(A, y, x0, xTrue, n, m, params, offsets, sizes, xNew, r);
                }
                break;
            case 3:
                cout << "\n=== Метод блочной релаксации ===\n";
                cout << "Введите размер блока: ";
                cin >> params.blockSize;
                params.maxIter = 1000000;
                cout << "Задайте диапазон omega через пробел(пример: 0.01 1.61): ";
                cin >> w0 >> w1;
                cout << "Задайте шаг omega (пример: 0.01): ";
                cin >> v;
                for (params.omega = w0; params.omega < w1; params.omega += v)
                {
                    BlockRelaxation(A, y, x0, xTrue, n, m, params, offsets, sizes, r);
                }
                break;
            case 0:
                
                exit = false;
                break;
            default:
                throw runtime_error("Неверный выбор метода. Введите 1 или 2.");
            }
        }
    }
    catch (exception& e)
    {
        cerr << "Ошибка: " << e.what() << endl;
    }
    
}