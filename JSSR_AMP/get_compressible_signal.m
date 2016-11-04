function [x supp_index partial_supp_index cardinal_par_supp_set ] = get_compressible_signal(N,K,sigma_1,sigma_0)

%Generates a compressible signal with K large cofficients 

x=zeros(1,N);
x(1:K)=sigma_1*randn(1,K) ;
x(K+1:N)=sigma_0*randn(1,N-K);
ind=randperm(N);
x=x(ind);
supp_index = find(ind<=K);
cardinal_par_supp_set = floor(K*.7);
partial_supp_index = supp_index(1:cardinal_par_supp_set);
x=x';

% stem(x);