#ifndef RBM_H
#define RBM_H

#include <complex>
#include <Eigen/Eigenvalues>
#include "./matrix.h"
#include "./lattice.h"
enum class neuralnet_id {RBM, FFN};
class Rbm{
public:
	Rbm() {};
	~Rbm() {};
	Rbm(const neuralnet_id& nn_id, const int& nsites, const int& start_pos_, const int& hidden_density)
  		{ init_nn( nn_id, nsites, start_pos_, hidden_density); }
	void init_nn(const neuralnet_id& nn_id, const int& nsites,const int& start_pos_, const int& hidden_density);
	void get_rbm_parameters(const RealVector& pvec);
	void get_vlayer(RealVector& sigma, const ivector& row) const;
	std::complex<double> get_rbm_amplitudes(const ivector& row) const;
	void get_derivatives(const RealVector& pvec, ComplexVector& derivatives, const ivector& row, const int& start_pos)const ;
	void compute_theta_table(const ivector& row)const;
	const int& num_vparams(void) const { return num_params_; }
	void update_theta_table(const int& spin, const int& tsite, const int& fsite) const;


private:
int num_vunits_;
int kernel_size_;
int num_params_;
int alpha;
int num_hunits_;
double param_b_;
int tot_sites_;
int start_pos;

Matrix kernel_;
Vector in_bias_;
Vector hl_bias_;
mutable Vector input_;
mutable Vector theta_;
};
#endif 
