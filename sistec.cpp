#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>

using namespace std;
using number = long double;
using matrixType = vector<vector<number>>;

class matrix
{
protected:
    matrixType data;
    int n, m;

public:
    matrix() : n(3), m(3)
    {
        data = eye(3);
    }

    matrix(int n, int m) : n(n), m(m)
    {
        data = zeros(n, m);
    }

    matrix(matrixType matrix_data) : n(matrix_data.size()), m(matrix_data[0].size())
    {
        data = matrix_data;
    }

    static matrixType zeros(int n, int m)
    {
        return matrixType(n, vector<number>(m, 0));
    }

    static matrixType ones(int n, int m)
    {
        return matrixType(n, vector<number>(m, 1));
    }

    static matrixType eye(int n)
    {
        matrixType id = matrixType(n, vector<number>(n, 0));
        for (int i = 0; i < n; i++)
        {
            id[i][i] = 1;
        }
        return id;
    }

    void fill(number value)
    {
        data = matrixType(n, vector<number>(m, value));
    }

    void print(const string &title)
    {
        cout << title << endl;
        for (const auto &row : data)
        {
            for (number value : row)
            {
                cout << fixed << setprecision(2) << value << '\t';
            }
            cout << endl;
        }
        cout << endl;
    }

    void input_fill()
    {
        cout << "Rellena la matriz " << n << 'x' << m << endl;
        for (auto &row : data)
        {
            for (number &cell : row)
            {
                float temp;
                cin >> temp;
                cell = static_cast<number>(temp);
            }
        }
    }

    // Operación por fila (cuenta desde 0 a n)
    void e(int destRow, int fromRow, number multiplier)
    {
        if (destRow >= m || fromRow >= n)
        {
            throw range_error("Error: Operación fuera de rango");
        }
        for (int i = 0; i < m; i++)
        {
            data[destRow][i] += data[fromRow][i] * multiplier;
        }
    }

    void set(int i, int j, number value)
    {
        data[i][j] = value;
    }

    number get(int i, int j)
    {
        return data[i][j];
    }

    vector<number> &operator[](int index)
    {
        return data[index];
    }

    // Modificar triangulando la matriz y retornar los factores
    matrix triangulated()
    {
        matrix A = *this;
        for (int i = 0; i < n - 1; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                number factor = A[j][i] / A[i][i];
                if (factor == 0)
                {
                    throw logic_error("Error: Hay un 0 en la diagonal principal");
                }
                for (int k = i; k < m; k++)
                {
                    A[j][k] -= factor * A[i][k];
                }
            }
            if (A[i + 1][i + 1] == 0.0)
            {
                throw logic_error("Error: La matriz es linealmente dependiente");
            }
        }
        return A;
    }
};

class AmpliateMatrix : public matrix
{
private:
    vector<number> b;
    matrixType L, U;

    // Resolver un sistema Ax = B por eliminación Gaussiana
    vector<number> solveXByGaussElim()
    {
        vector<number> X(n, 0);

        // Triangulación de la matriz
        matrix A = *this;
        for (int i = 0; i < n - 1; i++)
        {
            if (A[i][i] == 0)
            {
                throw logic_error("Error: Hay un 0 en la diagonal principal");
            }
            for (int j = i + 1; j < n; j++)
            {
                number factor = A[j][i] / A[i][i];
                for (int k = i; k < m; k++)
                {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }

        // Imprimir la matriz triangulada
        A.print("Matriz A triangulada: ");

        // Despeje hacia atrás
        for (int i = n - 1; i >= 0; i--)
        {
            number sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += A[i][j] * X[j];
            }
            if (A[i][i] == 0)
            {
                throw logic_error("Error: División por cero en despeje hacia atrás");
            }
            X[i] = (b[i] - sum) / A[i][i];
        }

        return X;
    }

    vector<number> solveXByGaussSeidel() {
        vector<number> X(n, 0);
        vector<number> nX(n, 0);
        number minError = 1e-6;  // Tolerancia para el error
        vector<number> errors(n, 9999);  // Inicializar errores grandes
        int maxIterations = 1000;  // Límite de iteraciones para evitar ciclos infinitos
        int iterationCount = 0;

        while (true) {
            bool converged = true;  // Bandera para verificar convergencia

            for (int i = 0; i < n; i++) {
                if (data[i][i] == 0) {
                    throw runtime_error("Error: División por cero detectada en la diagonal principal.");
                }

                number sum = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += data[i][j] * X[j];  // Utiliza los valores más recientes de X
                    }
                }

                nX[i] = (b[i] - sum) / data[i][i];  // Calcular nuevo valor de X[i]
                errors[i] = abs(nX[i] - X[i]);  // Calcular el error

                // cout << "Iteración " << iterationCount << ": X[" << i << "] anterior: " << X[i]
                 //    << ", nuevo: " << nX[i] << ", error: " << errors[i] << endl;

                if (errors[i] > minError) {
                    converged = false;
                }

                // Actualizar inmediatamente el valor de X[i] con el valor calculado
                X[i] = nX[i];
            }

            iterationCount++;

            // Si todos los errores están por debajo de minError o se alcanzó el límite de iteraciones
            if (converged || iterationCount >= maxIterations) {
                break;
            }
        }

        if (iterationCount >= maxIterations) {
            cout << "Advertencia: No se alcanzó la convergencia después de " << maxIterations << " iteraciones." << endl;
        }

        return X;
    }

    vector<number> solveXByLU(bool bChanged = true) {
        vector<number> X(n, 0);  // Vector de solución
        vector<number> Y(n, 0);  // Vector intermedio para L*Y = b

        if (bChanged) {
            // 1. Factorización LU de la matriz A (data)
            for (int i = 0; i < n; i++) {
                // Llenar la matriz U
                for (int j = i; j < n; j++) {
                    number sum = 0;
                    for (int k = 0; k < i; k++) {
                        sum += L[i][k] * U[k][j];
                    }
                    U[i][j] = data[i][j] - sum;
                }

                // Llenar la matriz L
                for (int j = i; j < n; j++) {
                    if (i == j) {
                        L[i][i] = 1;  // Diagonal de L es 1
                    } else {
                        number sum = 0;
                        for (int k = 0; k < i; k++) {
                            sum += L[j][k] * U[k][i];
                        }
                        L[j][i] = (data[j][i] - sum) / U[i][i];
                    }
                }
            }
        }
        

        // 2. Resolver L * Y = b (Sustitución hacia adelante)
        for (int i = 0; i < n; i++) {
            number sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * Y[j];
            }
            Y[i] = (b[i] - sum);
        }

        // 3. Resolver U * X = Y (Sustitución hacia atrás)
        for (int i = n - 1; i >= 0; i--) {
            number sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += U[i][j] * X[j];
            }
            X[i] = (Y[i] - sum) / U[i][i];
        }

        return X;
    }


public:
    enum SolveMethods
    {
        BY_GAUSS_ELIM,
        BY_GASUS_SEIDEL,
        BY_LU
    };
    SolveMethods solveMethod;
    AmpliateMatrix(matrixType matrix_data, vector<number> b) : matrix(matrix_data), b(b), solveMethod(SolveMethods::BY_LU), L(matrix::zeros(n, m)), U(matrix::zeros(n, m))
    {
    }

    void setB(vector<number> newB) {
        b = newB;
    }
    // Función para imprimir la matriz aumentada
    void print()
    {
        cout << "Matriz aumentada:" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                cout << fixed << setprecision(2) << data[i][j] << '\t';
            }
            cout << "| " << b[i] << endl;
        }
        cout << endl;
    }

    vector<number> solve(bool bChanged = true)
    {
        switch (solveMethod)
        {
        case SolveMethods::BY_GAUSS_ELIM:
            return solveXByGaussElim();
            break;
        case SolveMethods::BY_GASUS_SEIDEL:
            return solveXByGaussSeidel();
            break;
        case SolveMethods::BY_LU:
            return solveXByLU(bChanged);
            break; 
        default:
            break;
        }
        return vector<number>(0);
    }
};

// Sobrecarga del operador << para imprimir vectores
template <class T>
ostream &operator<<(ostream &os, const vector<T> &v)
{
    os << '(';
    for (int i = 0; i < v.size(); i++)
    {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << ')';
    return os;
}

int main()
{
    // Crear una matriz y un vector B
    matrixType A = {
        {4, -1, 0, 0},
        {-1, 4, -1, 0},
        {0, -1, 4, -1},
        {0, 0, -1, 3}
    };

    vector<number> b = {15, 10, 10, 10};


    AmpliateMatrix augMatrix(A, b);

    augMatrix.print();
    vector<number> X = augMatrix.solve(true);
    cout << "Solucion X: " << X << endl << endl;

    augMatrix.setB({12, 13, 14, 15});
    augMatrix.print();
    X = augMatrix.solve(true);
    cout << "Solucion X: " << X << endl;

    return 0;
}
