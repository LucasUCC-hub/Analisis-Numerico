#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <cmath>

using namespace std;
using matriz = vector<vector<double>>;

class Matriz {
protected:
    matriz datos;

public:
    int n, m;

    Matriz(int n, double relleno = 0) : n(n), m(n), datos(n, vector<double>(n, relleno)) {}
    Matriz(int n, int m, double relleno = 0) : n(n), m(m), datos(n, vector<double>(m, relleno)) {}
    Matriz(matriz A) : n(A.size()), m(A[0].size()), datos(A) {}
    Matriz(initializer_list<vector<double>> A) : n(A.size()), m(A.begin()->size()), datos(A) {}

    static matriz nula(int n, int m) { return matriz(n, vector<double>(m, 0)); }
    static matriz uno(int n, int m) { return matriz(n, vector<double>(m, 1)); }
    static matriz iden(int n) { matriz id(n, vector<double>(n, 0)); for (int i = 0; i < n; ++i) id[i][i] = 1; return id; }

    void rellenar(double value) { datos.assign(n, vector<double>(m, value)); }

    void imprimir(const string &titulo = "") const {
        if (!titulo.empty()) cout << titulo << endl;
        for (const auto &row : datos) {
            for (double value : row) {
                cout << fixed << setprecision(2) << value << '\t';
            }
            cout << endl;
        }
        cout << endl;
    }

    void rellenarPorTeclado() {
        cout << "Rellena la matriz " << n << 'x' << m << endl;
        for (auto &row : datos) {
            for (double &cell : row) {
                cin >> cell;
            }
        }
    }

    void e(int a, int de, double mult) {
        if (a >= n || de >= n) throw range_error("Error: Operación fuera de rango");
        for (int i = 0; i < m; ++i) datos[a][i] += datos[de][i] * mult;
    }

    double &sub(int i, int j) { return datos[i][j]; }
    const double &sub(int i, int j) const { return datos[i][j]; }

    Matriz triangulated() const {
        Matriz A = *this;
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double factor = A.sub(j, i) / A.sub(i, i);
                if (factor == 0) throw logic_error("Error: Hay un 0 en la diagonal principal");
                for (int k = i; k < m; ++k) A.sub(j, k) -= factor * A.sub(i, k);
            }
            if (A.sub(i + 1, i + 1) == 0) throw logic_error("Error: La matriz es linealmente dependiente");
        }
        return A;
    }

    Matriz punto(Matriz& A) {
        if (m != A.n) throw invalid_argument("Error: Dimensiones de matrices factor incompatibles");
        Matriz res(n, A.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < A.m; j++) {
                double sum = 0;
                for (int k = 0; k < A.n; k++) {
                    sum += sub(i, k) * A.sub(k, j); 
                }
                res.sub(i, j) = sum;
            }
        }
        return res;
    }


    Matriz inversa() {
        if (n != m) throw invalid_argument("Error: La matriz no es cuadrada");
        Matriz A = *this;
        Matriz I = Matriz::iden(n);

        // Aplicamos eliminación de Gauss-Jordan
        for (int i = 0; i < n; ++i) {
            // Verificamos que el pivote no sea cero
            if (A.sub(i, i) == 0) {
                bool swapped = false;
                for (int j = i + 1; j < n; ++j) {
                    if (A.sub(j, i) != 0) {
                        // Intercambiamos filas
                        swap(A.datos[i], A.datos[j]);
                        swap(I.datos[i], I.datos[j]);
                        swapped = true;
                        break;
                    }
                }
                // Alguna fila es linealmente dependiente
                if (!swapped) throw runtime_error("Error: La matriz es singular y no tiene inversa.");
            }

            // Normalizamos la fila para que el pivote sea 1
            double pivote = A.sub(i, i);
            for (int j = 0; j < n; ++j) {
                A.sub(i, j) /= pivote;
                I.sub(i, j) /= pivote;
            }

            // Hacemos ceros en las otras filas en la columna del pivote
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    double factor = -A.sub(j, i);
                    A.e(j, i, factor);
                    I.e(j, i, factor);
                }
            }
        }

        return I;
    }

    double normaFrobenius() {
        if (n != m) throw invalid_argument("Error: La matriz no es cuadrada");
        double sumaCuadrados = 0.0;

        // Sumar los cuadrados de cada elemento de la matriz
        for (const auto& fila : datos) {
            for (double valor : fila) {
                sumaCuadrados += valor * valor;
            }
        }

        // Retornar la raíz cuadrada de la suma de los cuadrados
        return sqrt(sumaCuadrados);
    }

    double normaInfinito() {
        double maxSuma = 0.0;

        // Iterar sobre cada fila de la matriz
        for (const auto& fila : datos) {
            double sumaFila = 0.0;

            // Sumar los valores absolutos de los elementos de la fila
            for (double valor : fila) {
                sumaFila += abs(valor);
            }

            // Actualizar el valor máximo de la suma
            if (sumaFila > maxSuma) {
                maxSuma = sumaFila;
            }
        }

        return maxSuma;
    }

    double normaUno() {
        double maxSuma = 0.0;

        // Iterar sobre cada columna
        for (int j = 0; j < m; ++j) {
            double sumaColumna = 0.0;

            // Sumar los valores absolutos de los elementos de la columna
            for (int i = 0; i < n; ++i) {
                sumaColumna += abs(datos[i][j]);
            }

            // Actualizar el valor máximo de la suma de columnas
            if (sumaColumna > maxSuma) {
                maxSuma = sumaColumna;
            }
        }

        return maxSuma;
    }

    double norma() {
        return normaInfinito();
    }

    int cond() {
        return norma() * inversa().norma();
    }

    double det() {
        if (n != m) throw invalid_argument("Error: La matriz no es cuadrada");

        Matriz A = *this;  // Hacer una copia de la matriz para no modificar la original
        double determinante = 1.0;

        for (int i = 0; i < n; ++i) {
            // Si el pivote es cero, intentar hacer un intercambio de filas
            if (A.sub(i, i) == 0) {
                bool filaIntercambiada = false;
                for (int j = i + 1; j < n; ++j) {
                    if (A.sub(j, i) != 0) {
                        swap(A.datos[i], A.datos[j]);
                        determinante *= -1;  // Intercambiar filas cambia el signo del determinante
                        filaIntercambiada = true;
                        break;
                    }
                }
                if (!filaIntercambiada) return 0.0;  // Si no hay fila no nula, el determinante es 0
            }

            // Escalar el pivote
            determinante *= A.sub(i, i);

            // Realizar eliminación de Gauss
            for (int j = i + 1; j < n; ++j) {
                double factor = A.sub(j, i) / A.sub(i, i);
                A.e(j, i, -factor);  // Utilizamos la función `e` para aplicar la operación por fila
            }
        }

        return determinante;
    }
};

class MatrizAumentada : public Matriz {
private:
    vector<double> b;
    matriz L, U;

    vector<double> resolverPorGaussElim() {
        Matriz A = *this;
        vector<double> X(n, 0);
        for (int i = 0; i < n - 1; ++i) {
            if (A.sub(i, i) == 0) throw logic_error("Error: Hay un 0 en la diagonal principal");
            for (int j = i + 1; j < n; ++j) {
                double factor = A.sub(j, i) / A.sub(i, i);
                for (int k = i; k < m; ++k) A.sub(j, k) -= factor * A.sub(i, k);
                b[j] -= factor * b[i];
            }

            cout << i + 1 << " PASO";
            A.imprimir("");
        }
        
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j) sum += A.sub(i, j) * X[j];
            if (A.sub(i, i) == 0) throw logic_error("Error: División por cero en despeje hacia atrás");
            X[i] = (b[i] - sum) / A.sub(i, i);
        }
        return X;
    }

    vector<double> resolverPorGaussSeidel() {
        vector<double> X(n, 0), nX(n, 0);
        const double minError = 1e-6;
        const int maxIterations = 1000;
        int iterationCount = 0;

        while (true) {
            bool converged = true;
            for (int i = 0; i < n; ++i) {
                if (sub(i, i) == 0) throw runtime_error("Error: División por cero detectada en la diagonal principal.");
                double sum = 0;
                for (int j = 0; j < n; ++j) if (j != i) sum += sub(i, j) * X[j];
                nX[i] = (b[i] - sum) / sub(i, i);
                if (abs(nX[i] - X[i]) > minError) converged = false;
                X[i] = nX[i];
            }
            if (converged || ++iterationCount >= maxIterations) break;
        }
        if (iterationCount >= maxIterations) cout << "Advertencia: No se alcanzó la convergencia después de " << maxIterations << " iteraciones." << endl;
        return X;
    }

    vector<double> resolverPorLU(bool bChanged = true) {
        vector<double> X(n, 0), Y(n, 0);
        if (bChanged) {
            L = Matriz::nula(n, n);
            U = Matriz::nula(n, n);
            for (int i = 0; i < n; ++i) {
                for (int j = i; j < n; ++j) {
                    double sum = 0;
                    for (int k = 0; k < i; ++k) sum += L[i][k] * U[k][j];
                    U[i][j] = sub(i, j) - sum;
                }
                for (int j = i; j < n; ++j) {
                    if (i == j) L[i][i] = 1;
                    else {
                        double sum = 0;
                        for (int k = 0; k < i; ++k) sum += L[j][k] * U[k][i];
                        L[j][i] = (sub(j, i) - sum) / U[i][i];
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < i; ++j) sum += L[i][j] * Y[j];
            Y[i] = b[i] - sum;
        }
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j) sum += U[i][j] * X[j];
            X[i] = (Y[i] - sum) / U[i][i];
        }
        return X;
    }

public:
    enum ResolverMetodos { GAUSS_ELIM, GAUSS_SEIDEL, LU };

    ResolverMetodos resolverMethod;

    MatrizAumentada(matriz A, vector<double> b, ResolverMetodos metodo = LU)
        : Matriz(A), b(move(b)), resolverMethod(metodo), L(Matriz::nula(n, n)), U(Matriz::nula(n, n)) {}

    void setB(vector<double> newB) { b = move(newB); }

    void imprimir() const {
        cout << "Matriz aumentada:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) cout << fixed << setprecision(2) << sub(i, j) << '\t';
            cout << "| " << b[i] << endl;
        }
        cout << endl;
    }

    vector<double> resolver(bool bChanged = true) {
        switch (resolverMethod) {
        case GAUSS_ELIM: return resolverPorGaussElim();
        case GAUSS_SEIDEL: return resolverPorGaussSeidel();
        case LU: return resolverPorLU(bChanged);
        default: throw invalid_argument("Método de solución no soportado.");
        }
    }
};


double numeroCondicion(const matriz &A) {
    double maxNorm = 0;
    for (const auto& row : A) {
        for (double val : row) {
            if (std::abs(val) > maxNorm) {
                maxNorm = std::abs(val);
            }
        }
    }

    return maxNorm;
}


template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
    os << '(';
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) os << ", ";
    }
    return os << ')';
}

int main() {
    Matriz A = {

        {2, 3, 4},
        {2, 1, 0},
        {8, 4, 4}
    };

    Matriz B = {

        {0, 6, 7},
        {1, 5, 7},
        {0, 2, 5}
    };

    A.imprimir();
    B.imprimir();
    A.punto(B).imprimir();

    int condicion = A.det();
    cout << "Condicion de A: " << condicion;

    return 0;
}

#include <iostream>

class Matrix {
private:
    int rows, cols;
    int* data;

public:
    // Constructor
    Matrix(int r, int c) : rows(r), cols(c) {
        data = new int[rows * cols]; // Crear un array unidimensional
        for (int i = 0; i < rows * cols; ++i) {
            data[i] = 0; // Inicializar en cero
        }
    }

    // Destructor
    ~Matrix() {
        delete[] data; // Liberar memoria
    }

    // Overloading () operator
    int& operator()(int x, int y) {
        return data[x * cols + y]; // Calcular el índice
    }

    const int& operator()(int x, int y) const {
        return data[x * cols + y];
    }

    // Print the matrix
    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                std::cout << data[i * cols + j] << " "; // Calcular el índice
            }
            std::cout << std::endl;
        }
    }
};

int main() {
    int r = 3, c = 3;
    Matrix m(r, c); // Crear una matriz de 3x3

    // Asignar valores a la matriz
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            m(i, j) = i * c + j; // Llenar la matriz con algunos valores
        }
    }

    // Imprimir la matriz
    m.print();

    // Acceso a un elemento específico
    std::cout << "Element at (1, 1): " << m(1, 1) << std::endl;

    return 0;
}
