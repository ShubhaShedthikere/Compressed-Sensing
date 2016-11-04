function A = get_encoding_matrix_bernoulli(m,n)

% Generates the encoding matrix whose entries are i.i.d random gaussian


A= sign(randn(m,n));
% 
% ind = find(A<0);
% A(ind)=0;


