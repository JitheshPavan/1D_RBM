// File: rbm.cpp
#include "rbm.h"
#include <algorithm> 
#include <vector>



void Rbm::init_nn(const neuralnet_id& nn_id, const int& nsites,const int& start_pos_, const int& hidden_density)
{

  num_vunits_ = 2*nsites;
  tot_sites_=nsites;
  alpha  =  hidden_density;
  num_hunits_  = alpha*num_vunits_;
  tot_sites_= nsites;
  kernel_.resize(num_vunits_,num_hunits_);
  in_bias_.resize(num_vunits_);
  hl_bias_.resize(num_hunits_);
  start_pos = start_pos_;
  num_params_ = in_bias_.size()+hl_bias_.size()+kernel_.size();//+num_sign_params_;  
  theta_.resize(num_vunits_,hidden_density);
}

void Rbm::get_rbm_parameters(const RealVector& pvec)
{
  int n = start_pos+in_bias_.size();
  int m = n+hl_bias_.size();
  int l = m+kernel_.size(); //in_bias_.size()+
  int p = l+1;
  for (int i=0; i<in_bias_.size(); ++i)
    *(in_bias_.data()+i) = pvec(start_pos+i);
  for (int i=0; i<hl_bias_.size(); i++)
    *(hl_bias_.data()+i) = pvec(n+i);
  for (int i=0; i<kernel_.size(); ++i)
    *(kernel_.data()+i) = pvec(m+i);
}


void Rbm::get_vlayer(RealVector& sigma, const ivector& row) const
{
  for (int i=0; i<tot_sites_; i++){
    if(row(i) == 1) sigma(i) = 1;
    else sigma(i) = 0;
    if(row(i+tot_sites_) == 1) sigma(i+tot_sites_) = 1;
    else sigma(i+tot_sites_) = 0;
  }  
}

void Rbm::compute_theta_table(const ivector& row) const
{
  RealVector sig;
  sig.resize(num_vunits_);
  get_vlayer(sig, row);
  theta_.setZero();
  int k = 0;  
  int counter = 0;
  for(int i=0; i<alpha; i++){
    for (int j=0; j<num_vunits_; j++){
      double b = hl_bias_(counter);
      double val = b+(sig.transpose()*kernel_.col(j+k));
      theta_(j,i)  = val;
      counter ++;
    }
    k += num_vunits_;
  }
}

std::complex<double> Rbm::get_rbm_amplitudes(const ivector& row) const
{
  RealVector sig;
  sig.resize(num_vunits_);
  get_vlayer(sig, row);
  double e=in_bias_.transpose()*sig;
  std::complex<double> F = {1.0,0.0};
  for (int i=0; i<alpha; i++) 
    for (int j=0; j<num_vunits_; j++)
      F = F*cosh(theta_(j,i));
  std::complex<double> sign = 1.0; 
  return exp(e)*sign*F; //exp(e)*F;
}
void Rbm::update_theta_table(const int& spin, const int& tsite, const int& fsite) const
{ 
  int ts,fs;

  if(spin==0){
    ts = tsite + tot_sites_;
    fs = fsite + tot_sites_;
  }
  else{
    ts = tsite;
    fs = fsite;
  }
  int k=0;
  for(int i=0; i<alpha; i++){
    for(int j=0; j<num_vunits_; j++){
      int d_ = k+j;
      double val = (kernel_(ts,d_) - kernel_(fs,d_));
      theta_(j,i) += val;
    }
    k += num_vunits_;
  }
}

void Rbm::get_derivatives(const RealVector& pvec, ComplexVector& derivatives, const ivector& row, const int& start_pos) const
{
  RealVector sigma;
  sigma.resize(num_vunits_);
  //get_visible_layer(sigma, row);
  get_vlayer(sigma, row);
  int n  = 0 +in_bias_.size();
  int m  = n+hl_bias_.size();
  int m1 = m+kernel_.size();
  int l  = m1+1;
  ComplexMatrix M3(in_bias_.size(), hl_bias_.size());
  ComplexVector a(hl_bias_.size());
  for (int i=0; i<in_bias_.size(); i++)
    derivatives(start_pos+i) = sigma(i);
  for (int i = 0; i < alpha; ++i){
    for (int j = 0; j < hl_bias_.size(); ++j){
      derivatives(n+j) = tanh(theta_(j,i));
      a(j) = derivatives(n+j);
    }
  }
  M3 = sigma*a.transpose();
  for (int i=0; i<kernel_.size(); i++)
    derivatives(m+i) = *(M3.data()+i);
  /*double sin_val = -PI*std::tan(PI*theta_sign_);
  derivatives(m1) = sin_val;
  for (int i = 0; i < sign_weights_.size(); ++i)
    derivatives(l+i) = sigma(i)*sin_val;*/
}


