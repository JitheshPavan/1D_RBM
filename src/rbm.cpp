// File: rbm.cpp
#include "rbm.h"
#include <algorithm> 
#include <vector>
#include <cassert>

void Rbm::init_nn( const int& nsites, const int& hidden_density)
{
  hblock_=nsites;
  num_sites_=nsites;
  alpha=2*hidden_density;
  num_hblocks_=2*hidden_density;//alpha

  num_hunits_= hblock_*num_hblocks_; 
  kernel_= Matrix::Zero(num_hunits_, 2*num_sites_);
  h_bias_= Vector::Zero(num_hunits_);

  num_kernel_params_= num_hblocks_*2*num_sites_;
  num_hbias_params_ = num_hblocks_;
  num_params_=num_kernel_params_+num_hbias_params_;

  tanh_.resize(num_hunits_);
  a_matrix.resize(num_sites_,num_sites_);
  a_matrix1.resize(num_sites_,num_sites_);
 /*
  std::cout << "Testing the code\n"<<std::endl;
  Vector pvec(num_params_);
  for(int i=0;i< num_params_;++i){
    if(i%3==0){
      pvec[i]=0.001;
    }
    if(i%3==1){
      pvec[i]=-0.05;
    }
    if(i%3==2){
      pvec[i]=0.02;
    }
  }
  ivector sigma(2*num_sites_);
  for(int i=0;i<2*num_sites_;++i){
    if(i%3==0 || i%3==1){
      sigma[i]=0;
    }
    else{
      sigma[i]=1;
    }
  }
  get_rbm_parameters(pvec);
  RealVector n2(1);
  compute_theta_table(sigma);
  std::cout << "kernel_"<< kernel_<<std::endl;
  std::cout << "h_bias_="<< h_bias_<<std::endl;

  std::cout << "theta"<< theta_<<std::endl;
  std::cout << "the amplitude="<<get_rbm_amplitudes(sigma);

  ComplexVector grad(num_params_);
  get_derivatives(n2,grad,sigma,0);

  std::cout<< grad.squaredNorm();

  std::cout << "gradients=\n";
  for(int i=0;i<num_params_;i++){
    std::cout << i<<":"<<grad(i)<<"\n"<<std::endl;
  }
  getchar();
  */
}

void Rbm::get_rbm_parameters(const RealVector& pvec)

{
  for(int k=0;k<num_hblocks_;k++){
    int temp1= k*num_sites_;
    for(int i=0;i< num_sites_;++i){
      for(int j=0;j<num_sites_;++j){
        int temp=i+j;
        while(temp>=num_sites_){
          temp=temp-num_sites_;
        }
        kernel_(i+temp1,j)=pvec[2*temp1+temp];
        kernel_(i+temp1,j+num_sites_)=pvec[2*temp1+num_sites_+temp];
      }
    }
  }
  for(int i=0;i<num_hblocks_;++i){
    for(int j=0;j<num_sites_;++j){
      h_bias_[j+(i*num_sites_)]=pvec[num_kernel_params_+i];
    }
  }
}


void Rbm::get_vlayer(RealVector& sigma, const ivector& row) const
{
  for (int i=0; i<num_sites_; i++){
    if(row(i) == 1) sigma(i) = 1;
    else sigma(i) = 0;
    if(row(i+num_sites_) == 1) sigma(i+num_sites_) = 1;
    else sigma(i+num_sites_) = 0;
  }  
}

void Rbm::compute_theta_table(const ivector& row) const
{
  RealVector sig;
  sig.resize(2*num_sites_);
  get_vlayer(sig, row);
  theta_.setZero();
//  theta_.size()=num_hunits_=num_hblock_*hblock_;
  theta_= kernel_ *sig+ h_bias_;
}

std::complex<double> Rbm::get_rbm_amplitudes(const ivector& row) const
{
  RealVector sig;
  sig.resize(2*num_sites_);
  get_vlayer(sig, row);
  std::complex<double> F = {1.0,0.0};
  for(int i=0;i<theta_.rows();++i){
    F *= cosh(theta_[i]);
  }
  return F; 
}
void Rbm::update_theta_table(const int& spin, const int& tsite, const int& fsite) const
{ 
  int ts,fs;

  if(spin==0){
    ts = tsite + num_sites_;
    fs = fsite + num_sites_;
  }
  else{
    ts = tsite;
    fs = fsite;
  }
  for(int i=0;i<num_hunits_;++i){
    theta_[i]+=kernel_(i,ts)-kernel_(i,fs);
  }
}

void Rbm::get_derivatives(const RealVector& pvec, ComplexVector& derivatives, const ivector& row, const int& start_pos) const
{
  RealVector sigma;
  sigma.resize(2*num_sites_);
  //get_visible_layer(sigma, row);
  get_vlayer(sigma, row);
  Vector tempvec;
  derivatives.setZero();
  //length of tanh_ = num_sites*num_hblocks_
  for(int i=0;i<num_hunits_;++i){
    tanh_[i]=tanh(theta_[i]);
  }

  for(int i=0;i<num_sites_;++i){
    for(int j=0;j<num_sites_;++j){
      int temp=i-j;
      while(temp <0){
        temp+= num_sites_;
      }
      a_matrix(i,j)=sigma[temp];
      a_matrix1(i,j)=sigma[temp+num_sites_];
    }
  }

  for(int i=0;i<num_hblocks_;++i){
    tempvec=tanh_.segment(i*num_sites_,num_sites_);
    der  =  a_matrix * tempvec;
    der1 = a_matrix1 * tempvec;
    for(int j=0;j<num_sites_;++j){
      if(std::abs(der(j))>1.0E-12) derivatives(i*2*num_sites_+j)=der(j);
      if(std::abs(der1(j))>1.0E-12) derivatives(i*2*num_sites_+num_sites_+j)=der1(j);
    }
  }

  for(int i=0;i<num_hblocks_;++i){
    derivatives(num_kernel_params_+i)= tanh_.segment(i*num_sites_,num_sites_).sum();
  }
}


