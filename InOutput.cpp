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
void InputDiagMatrix(Matrix& A, real*& xTrue, real*& x0, real*& y, int& n, int& m, string filename)
{

    try
    {
        ifstream fin(filename);

        if (!fin.is_open())
        {
            throw runtime_error("Не удалось открыть файл matrixDiag.txt");
        }

        fin >> n >> m;
        int sizes[9];
        CalculatingDimensionDiagonals(sizes, n, m);
        A.d4l = new real[sizes[0]]; for (int j = 0; j < sizes[0]; j++) fin >> A.d4l[j];
        A.d3l = new real[sizes[1]]; for (int j = 0; j < sizes[1]; j++) fin >> A.d3l[j];
        A.d2l = new real[sizes[2]]; for (int j = 0; j < sizes[2]; j++) fin >> A.d2l[j];
        A.d1l = new real[sizes[3]]; for (int j = 0; j < sizes[3]; j++) fin >> A.d1l[j];
        A.d0 = new real[sizes[4]]; for (int j = 0; j < sizes[4]; j++) fin >> A.d0[j];
        A.d1u = new real[sizes[5]]; for (int j = 0; j < sizes[5]; j++) fin >> A.d1u[j];
        A.d2u = new real[sizes[6]]; for (int j = 0; j < sizes[6]; j++) fin >> A.d2u[j];
        A.d3u = new real[sizes[7]]; for (int j = 0; j < sizes[7]; j++) fin >> A.d3u[j];
        A.d4u = new real[sizes[8]]; for (int j = 0; j < sizes[8]; j++) fin >> A.d4u[j];


        xTrue = new real[n];
        x0 = new real[n];
        y = new real[n];

        for (int i = 0; i < n; i++) fin >> xTrue[i];
        for (int i = 0; i < n; i++) fin >> x0[i];

        for (int i = 0; i < n; i++) fin >> y[i];

        fin.close();
    }
    catch (exception& e)
    {
        cout << "Ошибка: " << e.what() << endl;
        exit(0);
    }
}
void PrintMatrixStruct(const Matrix& A, int n, int m)
{
    cout << n << endl;
    int sizes[9];

    CalculatingDimensionDiagonals(sizes, n, m);

    for (int j = 0; j < sizes[0]; j++) cout << A.d4l[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[1]; j++) cout << A.d3l[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[2]; j++) cout << A.d2l[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[3]; j++) cout << A.d1l[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[4]; j++) cout << A.d0[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[5]; j++) cout << A.d1u[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[6]; j++) cout << A.d2u[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[7]; j++) cout << A.d3u[j] << ' ';
    cout << endl;

    for (int j = 0; j < sizes[8]; j++) cout << A.d4u[j] << ' ';
    cout << endl;
}

void PrintMatrix(real* A, int n)
{
    for (int j = 0; j < n; j++)
    {
        cout << setprecision(accuracy) << A[j] << ' ';
    }
    cout << endl;
}