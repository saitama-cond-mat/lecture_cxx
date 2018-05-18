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
    }
    return l * (x * legg(l,x) - legg(l-1,x))/(x*x-1);

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

int main() {

    //matrix_power();
    double x = 0.5;
    std::cout << legg(0, x) << " " << 1 << std::endl;
    std::cout << legg(1, x) << " " << x << std::endl;
    std::cout << legg(2, x) << " " << 0.5 * (3*x*x -1) << std::endl;

    std::cout << d_legg(0, x) << " " << 0 << std::endl;
    std::cout << d_legg(1, x) << " " << 1 << std::endl;
    std::cout << d_legg(2, x) << " " << 3 * x << std::endl;

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
