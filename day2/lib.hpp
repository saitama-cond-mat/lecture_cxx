#pragma once

#include <iostream>

inline int my_pow(int a, int n) {
    return 0;
}

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
