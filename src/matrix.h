#include <Eigen/Core>
#include<complex>

#ifndef MATRIX_H
#define MATRIX_H

const	double PI = 4.0*std::atan(1.0);
const	double TWO_PI = 2*PI;
const std::complex<double> II = std::complex<double>(0.0,1.0); 
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
using ivector = Eigen::VectorXi;
using RealVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;
using IntMatrix = Eigen::MatrixXi;
using RealMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ColVector = Eigen::VectorXcd;
using RowVector = Eigen::RowVectorXcd;
using double_Array = Eigen::ArrayXd;

#endif
