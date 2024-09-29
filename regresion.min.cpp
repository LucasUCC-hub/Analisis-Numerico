#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>

using namespace std;

template <class T, class A>
void map(T& v, const A& f) {
    transform(v.begin(), v.end(), v.begin(), f);
}

struct Punto {
    double x;
    double y;
    Punto(double x, double y) : x(x), y(y) {}
};

struct Dual {
    double a;
    double b;
    Dual(double a = 0, double b = 0) : a(a), b(b) {}
};

enum Func {
    LINEAL,
    EXPONENCIAL,
    POTENCIAL,
    CRECIMIENTO
};

Dual regresionLineal(vector<Punto>& v) {
    double n = v.size();
    double a0, a1, sumXY = 0, sumX = 0, sumY = 0, sumX2 = 0;

    for (Punto p : v) {
        sumX += p.x;
        sumY += p.y;
        sumXY += p.x * p.y;
        sumX2 += p.x * p.x;
    }

    a1 = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    a0 = (sumY / n) - a1 * (sumX / n);
    
    return Dual(a0, a1);
}

Dual predecir(vector<Punto>& v, Func modelo) {
    function<Punto(Punto)> inputTrans;
    function<Dual(Dual)> outputTrans;

    switch (modelo) {
    case EXPONENCIAL:
        inputTrans = [](Punto p){ return Punto(p.x, log(p.y)); };
        outputTrans = [](Dual d){ return Dual(exp(d.a), d.b); };
        break;
    case POTENCIAL:
        inputTrans = [](Punto p){ return Punto(log10(p.x), log10(p.y)); };
        outputTrans = [](Dual d){ return Dual(pow(10, d.a), d.b); };
        break;
    case CRECIMIENTO:
        inputTrans = [](Punto p){ return Punto(1.0 / p.x, 1.0 / p.y); };
        outputTrans = [](Dual d){ return Dual(1.0 / d.a, d.b / d.a); };
        break;
    case LINEAL:
    default:
        inputTrans = [](Punto p){ return p; };
        outputTrans = [](Dual d){ return d; };
        break;
    }

    // Aplicar la transformación de entrada
    map(v, inputTrans);

    // Obtener la regresión lineal
    Dual resultado = regresionLineal(v);

    // Aplicar la transformación de salida
    return outputTrans(resultado);
}

void imprimirModelo(const Dual& resultado, Func modelo, const string& function_name = "f") {
    switch (modelo) {
    case EXPONENCIAL:
        cout << function_name << "(x) = " << resultado.a << "e^" << resultado.b << "x" << endl;
        break;
    case POTENCIAL:
        cout << function_name << "(x) = " << resultado.a << "x^" << resultado.b << endl;
        break;
    case CRECIMIENTO:
        cout << function_name << "(x) = " << resultado.a << "x / (" << resultado.b << " + x)" << endl;
        break;
    case LINEAL:
    default:
        cout << function_name << "(x) = " << resultado.a << "x + " << resultado.b << endl;
        break;
    }
}

int main() {
    vector<Punto> muestra = {{2, 3}, {4, 5}, {7, 6}};

    // Predicción e impresión para regresión exponencial
    Dual resultadoExponencial = predecir(muestra, EXPONENCIAL);
    imprimirModelo(resultadoExponencial, EXPONENCIAL);

    // Predicción e impresión para regresión lineal
    Dual resultadoLineal = predecir(muestra, LINEAL);
    imprimirModelo(resultadoLineal, LINEAL, "g");

    return 0;
}
