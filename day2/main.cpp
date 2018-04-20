#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <complex>

class my_complex {
public:
    // Constructor
    my_complex(double re, double im) : _re(re), _im(im) {
        std::cout << " _re, _im = " << _re << " " << _im << std::endl;
    }

    double real() const {
        return _re;
    }

    double imag() const {
        return _im;
    }

    my_complex operator*(const my_complex& another_value) const {
        double re = this->_re * another_value.real() - this->_im * another_value.imag();
        double im = this->_re * another_value.imag() + this->_im * another_value.real();

        return my_complex(re, im);
    }

    my_complex& operator=(const my_complex& another_value) {
        _re = another_value.real();
        _im = another_value.imag();

        return *this;
    }

private:
    double _re, _im;
};


int main() {
    // Complex numbers
    std::complex<double> z(0.0, 1.0);
    std::complex<double> z2(1.0, 1.0);
    std::cout << z * z2 << std::endl;
    auto z3 = z * z2;
    std::cout << z3 << std::endl;

    std::cout << "Re " << z3.real() << std::endl;
    std::cout << "Im " << z3.imag() << std::endl;

    my_complex my_z(0.0, 1.0);
    my_complex my_z2(1.0, 1.0);
    std::cout << " my_z " << my_z.real() << " " << my_z.imag() << std::endl;

    my_complex my_z3 = my_z * my_z2;

    std::cout << " my_z3 " << my_z3.real() << " " << my_z3.imag() << std::endl;

    return 0;
}