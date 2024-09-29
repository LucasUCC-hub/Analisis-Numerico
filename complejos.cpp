#include <iostream>

using namespace std;

template <class T>
struct complex {
    T r;
    T i;

    complex(T a) : r(a), i(0) {}
    complex(T a, T b) : r(a), i(b) {}

    complex<T> operator+(const complex<T>& z) const {
        return complex<T>(r + z.r, i + z.i);
    }
};

ostream& operator<<(ostream& os, const complex<double>& z) {
    os << z.r << " + i" << z.i << endl;
    return os;
}

int main() {
    complex<double> c1(2, 3), c2(3, 4);
    cout << c1 + c2 << endl;
}