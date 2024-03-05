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
	Rbm( const int& nsites,const int& hidden_density)
  		{ init_nn( nsites, hidden_density); }
	void init_nn(const int& nsites,const int& hidden_density);
	void get_rbm_parameters(const RealVector& pvec);
	void get_vlayer(RealVector& sigma, const ivector& row) const;
	std::complex<double> get_rbm_amplitudes(const ivector& row) const;
	void get_derivatives(const RealVector& pvec, ComplexVector& derivatives, const ivector& row, const int& start_pos)const ;
	void compute_theta_table(const ivector& row)const;
	const int& num_vparams(void) const { return num_params_; }
	void update_theta_table(const int& spin, const int& tsite, const int& fsite) const;




private:
	int hblock_;
	int num_sites_;
	int num_hunits_;
	int num_hblocks_;
	int alpha;
	int num_kernel_params_;
	int num_hbias_params_;
	int num_params_;

	Matrix kernel_;
	Vector h_bias_;
	mutable Vector theta_;
	mutable Vector tanh_;
	mutable Matrix a_matrix;
	mutable Matrix a_matrix1;
	mutable Matrix der;
	mutable Matrix der1;
};	
#endif 
