#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <complex>
#include "./Eigen/Core"
#include "./Eigen/LU"
#include "./Eigen/Eigenvalues"

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

    matrix_power();

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
