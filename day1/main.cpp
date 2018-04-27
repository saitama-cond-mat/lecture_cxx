#include <iostream>
#include <vector>
#include <array>
#include <numeric>

void test_references() {
    // Reference
    int i = 0;
    std::cout << "i = " << i << std::endl;

    int& j = i;
    std::cout << "j = " << j << std::endl;

    j = 100;

    std::cout << "i = " << i << std::endl;
    std::cout << "j = " << j << std::endl;


    // Pointer (Do not use pointers in your practical code unless you know what you are actually doing...)
    int* p_i = &i;
    std::cout << "Make almost no sense!: the value of pointer = address in memory space = " << p_i << std::endl;
    std::cout << "One can see that the both point to the data at the same address " << &i << std::endl;
    std::cout << "The data stored at p_i is " << *p_i << std::endl;

    *p_i = 10000;
    std::cout << "The data has been overwritten !: i = " << i << std::endl;
}

void test_cast() {
    double d = 100.0;

    //More modern style since C++11: you can let the compiler infer the type of d.
    //auto d = 100.0;

    // You can use int i = ... (classic style before C++11)
    int i = static_cast<int>(d);
    // Again, the compiler can determine the type of i instead of you.
    // auto i = static_cast<int>(d);

    // no need for explicit conversion
    d = 99;
}

void test_array() {
    // since C++11 (may be statically allocated in the stack practically)
    std::array<double,4> a {10, 20, 30, 40};

    {
        int k = 0;
        {
            int k = 100;
            // k
        }
    }

    for (int i=0; i < a.size(); ++i) {
        // The index runs from 0 (not from 1).
        int k = 0;
        std::cout << " i = " << i << " " << a[i] << std::endl;
    }

    // Iterator
    for (auto it = a.begin(); it != a.end(); ++it) {
        std::cout << *it << std::endl;
    }

    // Since C++11 (range-based for statement)
    for (auto& v : a) {
        std::cout << "v = " << v << std::endl;
    }

    // This will crash your program (no safeguards)
    // a[1000] = 0.0;

    // Dynamic arrays (may be allocated in the heap)
    std::vector<double> b {10, 20, 30, 40};
    std::vector<double> c(4); // vector of length 4 (advanced: elements may not be initialized as you expect...)

    for (auto& v : b) {
        std::cout << v << std::endl;
    }
    std::cout << "Resizing b" << std::endl;
    b.resize(2);
    std::cout << "Size of b " << b.size() << std::endl;
    for (auto& v : b) {
        std::cout << v << std::endl;
    }

    // Many algorithms in the standard libraries
    b.resize(100);
    std::fill(b.begin(), b.end(), 10.0);
    // include numeric
    std::cout << " sum of b " << std::accumulate(b.begin(), b.end(), 0) << std::endl;
    c.resize(100);
    // Copying elements...
    std::copy(b.begin(), b.end(), c.begin());
}

int test_functions(int n, int& m) {
    std::cout << "in test_functions:  n = " << n << std::endl;
    std::cout << "in test_functions:  m = " << m << std::endl;

    // Write new data
    n = 1000;
    m = 1000;

    return m;
}

void test_const(const int& read_only_value) {
    std::cout << "read only value " << read_only_value;
}

void test_const2(const std::vector<double>& array) {
    // forbidden
    //array.resize(10);
    //array[0] = 1;
}


int main() {

    // Reference and pointer
    test_references();

    test_cast();

    test_array();

    // Scope of variables and functions
    {
        // The scope of n and m is between the two brackets.
        int n = 0, m = 0;
        auto r = test_functions(n, m);
        std::cout << "in main:  n = " << n << std::endl;
        std::cout << "in main:  m = " << m << std::endl;
        std::cout << "returned value = " << r << std::endl;
    }

    //test_const();

    //test_const2();


    return 0;
}