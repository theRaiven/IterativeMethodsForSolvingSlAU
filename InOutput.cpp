// InOutput.cpp
#include "Header.h"
void CalculatingDimensionDiagonals(int sizes[], int n, int m)
{

    int k{ 0 };
    for (int i = 3; i > 0; i--)
    {
        sizes[k++] = n - m - i - 1;
    }


    sizes[k++] = n - 1;
    sizes[k++] = n;
    sizes[k++] = n - 1;

    for (int i = 1; i <= 3; i++)
    {
        sizes[k++] = n - m - i - 1;
    }
}
void InputDiagMatrix(real**& A, real*& xTrue, real*& x0, real*& y, int& n, int& m, int& maxIter, string filename)
{

    try
    {
        ifstream fin(filename);

        if (!fin.is_open())
        {
            throw runtime_error("Не удалось открыть файл matrixDiag.txt");
        }

        fin >> n >> m;
        A = new real * [9];
        int sizes[9];
        CalculatingDimensionDiagonals(sizes, n, m);

        for (int i = 0; i < 9; i++)
        {
            A[i] = new real[sizes[i]];
            for (int j = 0; j < sizes[i]; j++)
            {
                fin >> A[i][j];
            }
        }

        xTrue = new real[n];
        x0 = new real[n];
        y = new real[n];

        for (int i = 0; i < n; i++) fin >> xTrue[i];
        for (int i = 0; i < n; i++) fin >> x0[i];
        fin >> maxIter;

        for (int i = 0; i < n; i++) fin >> y[i];

        fin.close();
    }
    catch (exception& e)
    {
        cout << "Ошибка: " << e.what() << endl;
        exit(0);
    }
}

