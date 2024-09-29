// POR LUCAS FABIANI 2024 - UCC Analisis Matematico
/*
    Este es un programa para resolver
    por diferentes metodos (segun convenga)
    raices y minimos de una funcion
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

#define MAX_IT 1000

using namespace std;

using func = function<double(double)>;

struct Intervalo
{
    double a, b;
    Intervalo(double a, double b) : a(a), b(b) {}
};

ostream &operator<<(ostream &os, const vector<double> &vec)
{
    os << '(';
    for (size_t i = 0; i < vec.size(); ++i)
    {
        os << vec[i];
        if (i != vec.size() - 1) os << ", ";
    }
    return os << ')';
}

class Funcion
{
public:
    func f;
    enum MetodosRaices { BISECCION, NEWTON_RAPSON, SECANTE, PUNTO_FIJO };
    Intervalo intBuscar;
    double minError;

    Funcion(func f, Intervalo I, double e) : f(f), intBuscar(I), minError(e) {}

    bool hayAsintota(const Intervalo &intervalo)
    {
        return abs(f(intervalo.b) - f(intervalo.a)) > 1e6 * abs(f(intervalo.a) - f(intervalo.b));
    }

    double derivadaEstimada(double x,  int orden = 1, double h = 1e-5)
    {
        auto df = [this, h](double x) -> double { return (f(x + h) - f(x)) / h; };
        if (orden == 1) return df(x);
        if (orden == 2) return (df(x + h) - df(x)) / h;
        return 0;
    }

    vector<Intervalo> obtenerIntervalosPorBiseccion(double paso = 0.5)
    {
        vector<Intervalo> raices;
        for (double x = intBuscar.a; x <= intBuscar.b; x += paso)
        {
            double xNuevo = x + paso;
            if (f(x) * f(xNuevo) < 0 || f(x) == 0) raices.emplace_back(x, xNuevo);
        }
        return raices;
    }

    vector<Intervalo> intervalosMinimos(double paso = 0.5) {
        vector<Intervalo> minimos;
        double xAnt, xPos;

        for (double x = intBuscar.a; x <= intBuscar.b; x += paso)
        {
            xAnt = x - paso;
            xPos = x + paso;
            if (f(x) < f(xAnt) && f(x) < f(xPos)) minimos.emplace_back(xAnt, xPos);
        }
        return minimos;
    }

    vector<double> raicesPorMetodo(function<double(double, double)> metodo)
    {
        vector<double> raices;
        for (const auto &intervalo : obtenerIntervalosPorBiseccion())
        {
            double a = intervalo.a, b = intervalo.b;
            for (int i = 0; i < MAX_IT; ++i)
            {
                double c = metodo(a, b);
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

    vector<double> raices(MetodosRaices metodo)
    {
        switch (metodo)
        {
        case BISECCION:
            return raicesPorMetodo([](double a, double b) { return (a + b) / 2; });
        case NEWTON_RAPSON:
            return raicesPorMetodo([this](double a, double b) {
                double x = (a + b) / 2;
                return x - f(x) / derivadaEstimada(x);
            });
        case SECANTE:
            return raicesPorMetodo([this](double a, double b) {
                return b - f(b) * (b - a) / (f(b) - f(a));
            });
        case PUNTO_FIJO:
            return raicesPorMetodo([this](double a, double b) {
                double raiz = (a + b) / 2;
                double l = 1 / (abs(derivadaEstimada(raiz)) + 1);
                return raiz - l * f(raiz);
            });
        default:
            throw runtime_error("Método no soportado.");
        }
    }

    vector<double> minimos() {
        vector<double> raices;
        auto metodo = [this](double a, double b) {
            double x = (a + b) / 2;
            return x - derivadaEstimada(x) / derivadaEstimada(x, 2);
        };

        for (const auto &intervalo : intervalosMinimos())
        {
            double a = intervalo.a, b = intervalo.b;
            for (int i = 0; i < MAX_IT; ++i)
            {
                double c = metodo(a, b);
                
                cout << "ITERACION " << i << endl;
                cout << "X: " << c << endl;
                cout << "e: " << abs(c - b) << endl;

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
};

// Main
int main()
{
    // auto f = [](double x) -> double { return 2 * x * x - 1; };
    // Funcion fun(f, Intervalo(-10, 10), 1e-6);

    // cout << setprecision(8);
    // cout << "RAICES" << endl;
    // vector<Funcion::MetodosRaices> metodos = {Funcion::BISECCION, Funcion::PUNTO_FIJO, Funcion::NEWTON_RAPSON, Funcion::SECANTE};
    // vector<string> nombres = {"Bisección", "Punto Fijo", "Newton-Raphson", "Secante"};

    // for (size_t i = 0; i < metodos.size(); ++i)
    // {
    //     cout << "Por " << nombres[i] << ": " << fun.raices(metodos[i]) << endl;
    // }

    auto f = [](double x) -> double { return 2 - ((x + 2) / exp(x)); };

    Funcion fun(f, Intervalo(-10, 10), 1e-3);

    cout << fun.minimos();

    return 0;
}
