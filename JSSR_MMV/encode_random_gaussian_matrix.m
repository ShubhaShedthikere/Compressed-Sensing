function A = encode_random_gaussian_matrix(m,n)

A_temp1 = randn(m,n);
A_temp = reshape(A_temp1,m*n,1);
small_index =  find(abs(A_temp)<1e-30);
A=A_temp;
A(small_index) = 1e-5;
A = reshape(A,m,n);

