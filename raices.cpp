#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <functional>
#include <iomanip>

#define MAX_IT 1000

using namespace std;
using func = function<double(double)>;

struct Intervalo
{
    double a;
    double b;
    Intervalo(double a, double b) : a(a), b(b) {}
};

// Sobrecarga del operador << para pair
ostream &operator<<(ostream &os, const Intervalo &I)
{
    os << '[' << I.a << ", " << I.b << ']';
    return os;
}

// Sobrecarga del operador << para vector
template <class T>
ostream &operator<<(ostream &os, const vector<T> &vec)
{
    os << '(';
    for (size_t i = 0; i < vec.size(); i++)
    {
        os << vec[i];
        if (i != vec.size() - 1)
        {
            os << ", ";
        }
    }
    os << ')';
    return os;
}

class Funcion
{
public:
    func f;
    enum MetodosRaices
    {
        BISECCION,
        NEWTON_RAPSON,
        SECANTE,
        PUNTO_FIJO
    };
    Intervalo intBuscar;
    double minError;

    Funcion(func f, Intervalo I, double e) : f(f), intBuscar(I), minError(e) {}

    bool hayAsintota(Intervalo &intervalo)
    {
        double delta = abs(f(intervalo.a) - f(intervalo.b));
        double fa = f(intervalo.a);
        double fb = f(intervalo.b);
        return abs(fb - fa) > 1e6 * delta;
    }

    double derivadaEstimada(double x, double h = 1e-5)
    {
        return (f(x + h) - f(x)) / h;
    }

    double derivadaEstimadaTVI(Intervalo &intervalo)
    {
        return (f(intervalo.b) - f(intervalo.a)) / (intervalo.b - intervalo.a);
    }

    // Buscar intervalos con raiz
    vector<Intervalo> obtenerIntervalosPorBiseccion(double paso = 0.5,
                                                    int prof = 10)
    {
        vector<Intervalo> raices;
        double x = intBuscar.a;
        for (double xNuevo = intBuscar.a + paso; xNuevo <= intBuscar.b; xNuevo += paso)
        {
            // Intersección con el eje x
            if (f(x) * f(xNuevo) < 0)
            {
                raices.push_back(Intervalo(x, xNuevo));
            }
            // El mismo punto es raiz
            if (f(x) == 0)
            {
                raices.push_back(Intervalo(xNuevo, xNuevo));
            }
            x = xNuevo;
        }

        return raices;
    }

    vector<double> raicesPorBiseccion()
    {
        vector<double> raices;
        vector<Intervalo> intervalos = obtenerIntervalosPorBiseccion();

        for (auto it = intervalos.begin(); it != intervalos.end(); it++)
        {
            double a = it->a;
            double b = it->b;
            double medio = (a + b) / 2;
            double error = abs(a - b);

            for (int i = 0; i < MAX_IT; ++i)
            {
                medio = (a + b) / 2;
                error = abs(a - b) / 2;
                if (error < minError)
                {
                    raices.push_back(medio);
                    break;
                }

                if (hayAsintota(*it))
                {
                    intervalos.erase(it);
                    cout << "Asintota en: " << *it << endl;
                    continue;
                }
                // La raiz está a la izquierda
                if (f(medio) * f(a) < 0)
                {
                    b = medio;
                }
                // La raiz está a la derecha
                else if (f(medio) * f(b) < 0)
                {
                    a = medio;
                }
                // La raiz está en el medio
                else if (f(medio) == 0)
                {
                    a = medio;
                    b = medio;
                }
            }
        }

        return raices;
    }

    vector<double> raicesPorNewtonRapson(func df)
    {
        vector<double> raices;
        vector<Intervalo> intervalos = obtenerIntervalosPorBiseccion();
        for (auto &intervalo : intervalos)
        {
            double x = (intervalo.a + intervalo.b) / 2;
            for (int i = 0; i < MAX_IT; ++i)
            {
                double xNuevo = x - (f(x) / df(x));
                if (abs(xNuevo - x) < minError)
                {
                    raices.push_back(x);
                    break;
                }
                x = xNuevo;
            }
        }
        return raices;
    }

    vector<double> raicesPorSecante()
    {
        vector<double> raices;
        vector<Intervalo> intervalos = obtenerIntervalosPorBiseccion();
        for (auto &intervalo : intervalos)
        {
            double a = intervalo.a;
            double b = intervalo.b;
            for (int i = 0; i < MAX_IT; ++i)
            {
                double c = b - (f(b) * (b - a)) / (f(b) - f(a));
                if (abs(c - b) < minError)
                {
                    raices.push_back(c);
                    break;
                }
                a = b;
                b = c;
            }
        }
        return raices;
    }

    vector<double> raicesPorPuntoFijo()
    {
        vector<double> raices;
        vector<Intervalo> intervalos = obtenerIntervalosPorBiseccion();

        for (auto &intervalo : intervalos)
        {
            double raiz = (intervalo.a + intervalo.b) / 2;
            double l = 1 / (abs(derivadaEstimada(raiz)) + 1); // Parámetro de relajación
            auto g = [this, l, raiz](double x)
            { return x - l * f(x); };

            for (size_t i = 0; i < MAX_IT; i++)
            {
                double nuevaRaiz = g(raiz);

                if (abs(nuevaRaiz - raiz) < minError)
                {
                    raices.push_back(nuevaRaiz);
                    break;
                }

                raiz = nuevaRaiz;
            }
        }

        return raices;
    }

    vector<double> raices(MetodosRaices metodo = MetodosRaices::NEWTON_RAPSON)
    {
        switch (metodo)
        {
        case BISECCION:
            return raicesPorBiseccion();
        case NEWTON_RAPSON:
            return raicesPorNewtonRapson([this](double x) { return derivadaEstimada(x); });
        case PUNTO_FIJO:
            return raicesPorPuntoFijo();
        case SECANTE:
            return raicesPorSecante();
        default:
            throw runtime_error("Método no soportado.");
        }
    }
};

// Main
int main()
{
    auto f = [](double x) -> double
    { return 2 * x * x - 1; };

    Funcion fun(f, Intervalo(-10, 10), 1e-6);

    cout << setprecision(8);
    cout << "RAICES" << endl;
    cout << "Por Biseccion: " << fun.raices(Funcion::BISECCION) << endl;
    cout << "Por Punto Fijo: " << fun.raices(Funcion::PUNTO_FIJO) << endl;
    cout << "Por Newton-Raphson: " << fun.raices(Funcion::NEWTON_RAPSON) << endl;
    cout << "Por Secante: " << fun.raices(Funcion::SECANTE) << endl;

    return 0;
}
