#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <complex>
#include "./Eigen/Core"
#include "./Eigen/LU"
#include "./Eigen/Eigenvalues"

#include <utility>

double
legg_next(int l, const double &x, const double &Pl, const double &Plm1) {
    return ((2 * l + 1) * x * Pl - l * Plm1) / (l + 1);
}

double
legg(int l, double x) {
    double p0(1);
    double p1(x);

    if (l == 0) {
        return p0;
    }

    int n = 1;

    while (n < l) {
        std::swap(p0, p1);
        p1 = legg_next(n, x, p0, p1);
        ++n;
    }
    return p1;
}

double
d_legg(int l, double x) {
    if (l == 0) {
        return 0.0;
    } else if (l == 1) {
        return 1.0;
    }

    return ((2*l-1) * (legg(l-1,x) + x*d_legg(l-1,x)) - (l-1) * d_legg(l-2,x))/l;
    //return l * (x * legg(l,x) - legg(l-1,x))/(x*x-1);
}

double bisect(std::function<double(double)>& f, double x0, double x1, double aps = 1e-12) {
    double x_min = x0;
    double x_max = x1;

    assert(f(x_max) * f(x_min) < 0);

    while(x_max - x_min > aps) {
        double x_mid = 0.5*(x_max + x_min);

        if (f(x_mid) * f(x_min) > 0) {
            x_min = x_mid;
        } else {
            x_max = x_mid;
        }
    }

    return 0.5*(x_max + x_min);
}

std::vector<double> find_zeros(std::function<double(double)>& f) {
    //make double exponential mesh
    int N = 10000;
    std::vector<double> xs(N);
    double t_max = 3;
    double dt = 2*t_max/(N-1);
    for(int i=0; i<N; ++i) {
        xs[i]= std::tanh(0.5*M_PI*std::sinh(dt * i - t_max));
    }

    std::vector<double> zeros;
    for (int i=0; i<N-1; ++i) {
        // a sign change found
        if (f(xs[i]) * f(xs[i+1]) < 0)  {
            zeros.push_back(bisect(f, xs[i], xs[i+1] ));
        }
    }

    return zeros;
}

void matrix_power() {
    int P = 4;
    int N = 2;

    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(N,N);
    Eigen::MatrixXd rnd = (tmp + tmp.transpose()).eval();

    Eigen::EigenSolver<Eigen::MatrixXd> solver(rnd);

    Eigen::MatrixXcd D = Eigen::MatrixXcd::Zero(N,N);
    for (int i=0; i < N; ++i) {
        D(i,i) = std::pow(solver.eigenvalues()[i], P);
    }

    Eigen::MatrixXcd U = solver.eigenvectors();
    Eigen::MatrixXd r = (U * D * U.transpose()).real();

    std::cout << r << std::endl;
    std::cout << rnd * rnd * rnd * rnd << std::endl;
}

inline
void
assert_near(double x, double x2, double tol = 1e-8) {
    assert(std::abs(x-x2) < tol);
}

int main() {

    double x = 0.5;
    assert_near(legg(0, x), 1);
    assert_near(legg(1, x), x);
    assert_near(legg(2, x), 0.5 * (3*x*x -1));


    auto diff = [](int l, double x) {
        return (x * x - 1) * d_legg(l, x) - l * (x * legg(l, x) - legg(l-1, x));
    };

    for (auto x : std::vector<double>{-1.0, -0.1, 0.0, 0.1, 1.0}) {
        for (int l = 1; l < 10; ++l) {
            assert_near(diff(l,x), 0.0);
        }
    }
    //std::cout << legg(1, x) << " " << x << std::endl;
    //std::cout << legg(2, x) << " " << 0.5 * (3*x*x -1) << std::endl;

    //std::cout << d_legg(0, x) << " " << 0 << std::endl;
    //std::cout << d_legg(1, x) << " " << 1 << std::endl;
    //std::cout << d_legg(2, x) << " " << 3 * x << std::endl;
    //for (int l=0; l<10; ++l) {
        //std::cout << d_legg(l, 1.0) << " " << 0.5 * l * (l+1) << std::endl;
    //}

    int l = 10;
    std::function<double(double)> f = [&](double x){return legg(l, x);};
    auto zeros = find_zeros(f);

    for (auto x : zeros) {
        std::cout << x << " " << legg(l, x) << std::endl;
    }
    //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M;
    /*
    Eigen::MatrixXd M(100, 100);
    std::cout << M.rows() << " " << M.cols() << std::endl;
    M.resize(2, 2);
    std::cout << M.rows() << " " << M.cols() << std::endl;
    M.setZero();
    M(0,0) = 100.0;
    M(1,1) = 100.0;

    Eigen::MatrixXd inv_M = M.inverse();
    std::cout << inv_M << std::endl;
    std::cout << inv_M * M << std::endl;
     */


    //std::cout << solver.eigenvectors().conjugate().transpose() << " : " << D << " : " << solver.eigenvectors() << std::endl;

    //std::cout << rnd * rnd << std::endl;

    //Eigen::EigenSolver<Eigen::MatrixXd> solver(M;)

    //for (auto j = 0; j < 2; ++j) {
        //for (auto i = 0; i < 2; ++i) {
            //std::cout << i << " " << j << " " << M(i,j) << std::endl;
        //}
    //}

    //std::cout << M << std::endl;
    //assert(M.rows() >= 2);

    //auto M2 = M;

    //auto M3 = (M * M2).eval();
    //auto M3 = static_cast<Eigen::MatrixXd>(M + 2 * M2);

    //Eigen::MatrixXd M3;
    //M3 = M * M2;
    //auto M4 = M * M2;

    //std::cout << M3;
    //Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(N,N);
    //std::cout << (tmp + tmp.transpose()).eval() << std::endl;
}
