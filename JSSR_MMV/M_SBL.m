%function [xrec] = M_SBL(y,A,m,n,sigma_Z,l)
 clear all;
 clc;
 [x,y,n,k,m,l,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,support_set,w] = get_input_data_MMV();
num_iter_count=100;
sigma_z = sum(sigma_Z)/l;

var_z =sigma_z*sigma_z;
g = ones(n,1)*1;

%compute the posterior moments
for iter_count =1:num_iter_count
Gamma = diag(g);
sigma_t =  var_z*eye(m) +A*Gamma*A';

inverse_sigma_t = inv(sigma_t);
cov= Gamma - Gamma*A'*inverse_sigma_t*A*Gamma;

M = Gamma*A'*inverse_sigma_t*y;

%find new g

mu_temp = M.*M;
diag_cov_elements = cov.*eye(n);
inv_g = 1./g;
g_old=g;
if(l>1)
g=(sum(mu_temp')/l)'+ sum(diag_cov_elements)';
else
g = mu_temp + sum(diag_cov_elements)'
end

%var_z= (norm(y-A*M,'fro')^2/l)/(m-n + sum(diag_cov_elements)*(1./g));
xrec(:,iter_count) = M;
%mse(iter_count) =norm(x-M,'fro')^2 /(norm(M,2)^2);
end


