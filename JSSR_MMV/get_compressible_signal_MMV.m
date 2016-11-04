function [x support_set ] = get_compressible_signal_MMV(N,K,sigma_1,sigma_0,L)

%Generates a compressible signal with K large cofficients 

x=zeros(N,L);
x((1:K),:)=sigma_1*randn(K,L) ;
x((K+1:N),:)=sigma_0*randn(N-K,L);
ind=randperm(N);

for i=1:L
    x(:,i)= x(ind,i);
end

supp_index = find(ind<=K);
support_set=zeros(1,N);
support_set(supp_index)=1;



% stem(x);