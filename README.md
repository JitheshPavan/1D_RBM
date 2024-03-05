# RBM 
->Variational Monte Carlo with Stochastic Reconfiguration and RBM as NQS.</br>
RBM output is given as </br>
 W*x +b where w is the kernel and b is bias.</br>
RBM kernel matrix dimension is : 2 &times alpha*n_sites*,2*n_sites</br>
a single block of RBM is goes as : </br>
 [w1,w2,w3.... wn,  w(n+1),.....w(2n)] </br>
 [w2,w3,w4.....w1,  w(n+2),.....w(n+1)] </br>
 .
 .
 .
 [wn,w1,w2,....w(n-1), w(2n),w(n+1),.....w(2n-1)]
 (These blocks are repeated 2*alpha times for different sets of parameters).</br>
 We add a bias value for every block of the matrix. 
 Here n is num_sites_. This implements spin flip symmetry and translational symmetry.
the no of bias parameters = 2*alpha. That is, we add the same number for every n_sites output. 
 
