#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <complex>

#include "lib.hpp"


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