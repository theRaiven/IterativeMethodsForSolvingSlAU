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

        if (!(fin >> n >> m))
            throw runtime_error("Ошибка чтения размеров матрицы");

        if (n <= 0 || m < 0)
            throw runtime_error("Некорректные размеры n или m");

        int sizes[9];
        CalculatingDimensionDiagonals(sizes, n, m);
        auto readArray = [&](real*& arr, int size, const string& name)
        {
            arr = new real[size];
            for (int i = 0; i < size; i++)
            {
                if (!(fin >> arr[i]))
                    throw runtime_error("Недостаточно данных для " + name);
            }
        };

        readArray(A.d4l, sizes[0], "d4l");
        readArray(A.d3l, sizes[1], "d3l");
        readArray(A.d2l, sizes[2], "d2l");
        readArray(A.d1l, sizes[3], "d1l");
        readArray(A.d0, sizes[4], "d0");
        readArray(A.d1u, sizes[5], "d1u");
        readArray(A.d2u, sizes[6], "d2u");
        readArray(A.d3u, sizes[7], "d3u");
        readArray(A.d4u, sizes[8], "d4u");

        for (int i = 0; i < sizes[4]; i++)
        {
            if (A.d0[i] == 0)
            {
                throw runtime_error("Система вырождена: элемент d0[" + to_string(i) + "] = 0");
            }
        }


        xTrue = new real[n];
        x0 = new real[n];
        y = new real[n];

        auto readVector = [&](real* vec, int size, const string& name)
        {
            for (int i = 0; i < size; i++)
            {
                if (!(fin >> vec[i]))
                    throw runtime_error("Недостаточно данных для вектора " + name);
            }
        };

        readVector(xTrue, n, "xTrue");
        readVector(x0, n, "x0");
        readVector(y, n, "y");

        fin.close();
    }
    catch (exception& e)
    {
        cout << "Ошибка: " << e.what() << endl;
        exit(0);
    }
}
void WriteResultsToFile(const std::string& methodName, real w, int iter, real normR,
    real* x, real* xTrue, int n, const std::string& filename)
{
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open())
    {
        std::cerr << "Ошибка открытия файла " << filename << std::endl;
        return;
    }

    out << "=== " << methodName << " ===\n";
    out << "Достигнуто ";
    if (normR < 1e-15)
        out << "требуемое точность решения.\n";
    else
        out << "максимальное число итераций.\n";

    out << "w = " << std::setprecision(15) << w
        << ", Iter = " << iter
        << ", normR = " << std::setprecision(15) << normR << "\n";

    out << "i        x_i\n";
    for (int i = 0; i < n; i++)
        out << std::setprecision(15) << x[i] << "\n";

    if (xTrue != nullptr)
    {
        out << "i        x* _i - x_i\n";
        for (int i = 0; i < n; i++)
            out << std::setprecision(15) << xTrue[i] - x[i] << "\n";
    }

    out << "\nРезультат записан в файл: " << filename << "\n\n";
    out.close();
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