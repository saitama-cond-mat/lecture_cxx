#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <complex>

template<typename T>
class my_complex {
public:
    // Constructor
    my_complex(const my_complex& other) {
        std::cout << "Copy constructor " << std::endl;
        _re = other.real();
        _im = other.imag();
    }

    my_complex(T re, T im) : _re(re), _im(im) {
        std::cout << " _re, _im = " << _re << " " << _im << std::endl;
    }

    T real() const {
        return _re;
    }

    T imag() const {
        return _im;
    }

    my_complex operator*(const my_complex& another_value) const {
        T re = this->_re * another_value.real() - this->_im * another_value.imag();
        T im = this->_re * another_value.imag() + this->_im * another_value.real();

        return my_complex(re, im);
    }

    my_complex& operator=(const my_complex& another_value) {
        _re = another_value.real();
        _im = another_value.imag();

        std::cout << "operator=" << std::endl;

        return *this;
    }

private:
    T _re, _im;
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, const my_complex<T>& value) {
    stream << "( " << value.real() << " , " << value.imag() << " )";
    return stream;
}

template<typename T, typename T2>
T my_function(const T& a, const T2& b) {
    return a * b;
}


int main() {
    // Complex numbers
    std::complex<double> z(0.0, 1.0);
    std::complex<double> z2(1.0, 1.0);
    std::cout << z * z2 << std::endl;
    auto z3 = z * z2;
    std::cout << z3 << std::endl;

    std::cout << "Re " << z.real() << std::endl;
    std::cout << "Im " << z.imag() << std::endl;
    std::cout << "Re " << z3.real() << std::endl;
    std::cout << "Im " << z3.imag() << std::endl;

    my_complex<double> my_z(0.0, 1.0);
    my_complex<double> my_z2(1.0, 1.0);
    std::cout << " my_z " << my_z.real() << " " << my_z.imag() << std::endl;

    my_complex<double> my_z3 = my_z * my_z2;

    std::cout << " my_z3 " << my_z3.real() << " " << my_z3.imag() << std::endl;

    (std::cout << my_z3) << std::endl;

    //my_complex empty(my_z3);
    my_complex<double> empty2 = my_z3;

    std::cout << "my_func " << my_function<double,int>(1.0, 2) << std::endl;

    return 0;
}