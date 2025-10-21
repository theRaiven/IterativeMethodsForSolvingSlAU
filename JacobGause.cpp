#include "Header.h"
#include <cmath>
using namespace std;
void IterativeMethod(real** A, real* f, real* x, real* xTrue, int n, int m, real omega, int maxIter, bool isJacobi)
{
    if (!A || !f || !x || !xTrue)
    {
        cerr << "Ошибка: один из входных массивов равен nullptr!" << endl;
        return;
    }
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
    
    int sizes[9];
    CalculatingDimensionDiagonals(sizes, n, m);

    int offsets[]{ -4 - m, -3 - m, -2 - m, -1, 0, 1, 2 + m, 3 + m, 4 + m };

    real normR{ 0.0 };
    int iter{ 0 };
    for (iter = 1; iter < maxIter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            real sum{ 0.0 };

            for (int j = 0; j < n; j++)
            {
                int diff{ j - i };
                int d { -1 };

                if (diff <= -2 - m)
                {
                    d = diff + 4 + m;
                }
                else if (diff >= -1 && diff <= 1)
                {
                    d = diff + 4;
                }
                else if (diff >= 2 + m)
                {
                    d = diff - 2 - m + 6;
                }
                if (d < 0 || d > 8) continue;

                int offset{ 0 };
                if (d <= 2) offset = -4 - m + d;      // левая группа
                else if (d <= 5) offset = d - 4;      // центральная группа
                else offset = d - 6 + 2 + m;          // правая группа

                int i_min{ max(0, -offset) };
                int idx{ i - i_min };

                if (idx < 0 || idx >= sizes[d]) continue; // элемент вне диапазона
                real aij{ A[d][idx] };

                if (isJacobi)
                {
                    sum += aij * x[j];
                }
                else
                {
                    sum += aij * (j < i ? x[j] : x[j]);
                }
            }
            
            real x_next = (f[i] - sum) / A[4][i];
            if (isnan(x_next) || isinf(x_next))
            {
                cerr << "Ошибка: вычислен NaN/Inf на итерации " << iter << " в элементе " << i << endl;
                delete[] xNew;
                delete[] r;
                return;
            }

            
            if (isJacobi)
                xNew[i] = x[i] + omega * (x_next);
            else
                x[i] = x[i] + omega * (x_next);
        }

        if (isJacobi)
        {
            for (int i = 0; i < n; i++)
                x[i] = xNew[i];
        }
        // === вычисление относительной невязки ===
        normR = 0.0;
        real normF = 0.0;

        for (int i = 0; i < n; i++)
        {
            real Ax{ 0.0 };

            for (int j = 0; j < n; j++)
            {
                int diff{ j - i };
                int d{ -1 };

                if (diff <= -2 - m)
                    d = diff + 4 + m;          // левая группа
                else if (diff >= -1 && diff <= 1)
                    d = diff + 4;              // центральная группа
                else if (diff >= 2 + m)
                    d = diff - 2 - m + 6;      // правая группа

                if (d < 0 || d > 8) continue;

                int offset = offsets[d];
                int i_min = max(0, -offset);
                int idx = i - i_min;
                if (idx < 0 || idx >= sizes[d]) continue; // элемент вне диапазона

                Ax += A[d][idx] * x[j];
            }

            r[i] = f[i] - Ax;
            normR += r[i] * r[i];
            normF += f[i] * f[i];
        }


        normR = sqrt(normR / normF);

        if (normR < relativeEPS<real>())
        {
            cout << "Достигнута требуемая точность решения." << endl;
            break;
        }

        if (iter >= maxIter-1)
        {
            cout << "Достигнуто максимальное число итераций." << endl;
        }

    }

    cout << std::scientific << std::setprecision(15);
    cout << "w = " << omega << ", Iter = " << iter << ", normR = " << normR << endl;
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

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
    }
    delete[] xNew;
    delete[] r;
}
void BlockRelaxation(real** A, real* f, real* x, real* xTrue, int n, int m, real omega, int maxIter, int blockSize)
{
    if (!A || !f || !x || !xTrue)
    {
        cerr << "Ошибка: один из входных массивов равен nullptr!" << endl;
        return;
    }
    if (n <= 0 || blockSize <= 0 || blockSize > n)
    {
        cerr << "Ошибка: некорректный размер блока!" << endl;
        return;
    }

    real* r = new real[n];
    real normR{ 0.0 };

    int sizes[9];
    CalculatingDimensionDiagonals(sizes, n, m);
    int offsets[]{ -4 - m, -3 - m, -2 - m, -1, 0, 1, 2 + m, 3 + m, 4 + m };

    int iter = 0;
    for (iter = 1; iter <= maxIter; iter++)
    {
        for (int i0 = 0; i0 < n; i0 += blockSize)
        {
            int i1{ min(i0 + blockSize, n) };
            int bSize{ i1 - i0 };

            for (int i = i0; i < i1; i++)
            {
                real sum{ 0.0 };
                for (int j = 0; j < n; j++)
                {
                    int diff{ j - i };
                    int d{ -1 };
                    if (diff <= -2 - m)
                        d = diff + 4 + m;
                    else if (diff >= -1 && diff <= 1)
                        d = diff + 4;
                    else if (diff >= 2 + m)
                        d = diff - 2 - m + 6;

                    if (d < 0 || d > 8) continue;

                    int offset{ offsets[d] };
                    int i_min{ max(0, -offset) };
                    int idx{ i - i_min };
                    if (idx < 0 || idx >= sizes[d]) continue;

                    real aij{ A[d][idx] };
                    sum += aij * x[j];
                }
                r[i] = f[i] - sum;
            }

            // A_ii * Δx_i = ω * r_i
            for (int i = i0; i < i1; i++)
            {
                real delta = (omega * r[i]) / A[4][i];
                x[i] += delta;
            }
        }

        normR = 0.0;
        real normF = 0.0;

        for (int i = 0; i < n; i++)
        {
            real Ax{ 0.0 };
            for (int j = 0; j < n; j++)
            {
                int diff{ j - i };
                int d{ -1 };
                if (diff <= -2 - m)
                    d = diff + 4 + m;
                else if (diff >= -1 && diff <= 1)
                    d = diff + 4;
                else if (diff >= 2 + m)
                    d = diff - 2 - m + 6;

                if (d < 0 || d > 8) continue;

                int offset{ offsets[d] };
                int i_min{ max(0, -offset) };
                int idx{ i - i_min };
                if (idx < 0 || idx >= sizes[d]) continue;

                Ax += A[d][idx] * x[j];
            }

            real ri{ f[i] - Ax };
            normR += ri * ri;
            normF += f[i] * f[i];
        }

        normR = sqrt(normR / normF);
        

        if (normR < relativeEPS<real>())
        {
            cout << "Достигнута требуемая точность решения." << endl;
            break;
        }
    }
    cout << fixed << setprecision(2) << "w = " << omega<< ", Iter = " << iter;
    cout << std::scientific << std::setprecision(15) << ", normR = " << normR << endl;
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

    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
    }
    delete[] r;
}


void WorkingWithIterMethods()
{
	int n, m, maxIter;
    real** A = nullptr;
    real* xTrue = nullptr;
    real* x0 = nullptr;
    real* y = nullptr;

    try
    {
        InputDiagMatrix(A, xTrue, x0, y, n, m, maxIter, "Matrix.txt");
        cout << endl << n << ' ' << m << ' ' << maxIter << endl;
        PrintMatrix(A, n, m);
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
                IterativeMethod(A, y, x0, xTrue, n, m, w, maxIter, true);
            }
            break;

        case 2:
            cout << "\n=== Метод Гаусса–Зейделя ===\n";
            for (float w = 0.01; w < 1.70; w += 0.1)
            {
                IterativeMethod(A, y, x0, xTrue, n, m, w, maxIter, false);
            }
            break;
        case 3:
            cout << "\n=== Метод блочной релаксации ===\n";
            int blockSize;
            cout << "Введите размер блока: ";
            cin >> blockSize;
            for (float w = 0.01; w < 2; w += 0.1)
            {
                BlockRelaxation(A, y, x0, xTrue, n, m, w, maxIter, blockSize);
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