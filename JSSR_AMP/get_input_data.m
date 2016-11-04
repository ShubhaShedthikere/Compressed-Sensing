function [x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter, num_outer_iter,support_set] = get_input_data()

rho = 0.33;
delta = 0.3

n=5000;
m=floor(delta*n)
k=floor(rho*m);
s=k/n;
num_inner_iter=50;
num_outer_iter=30;

                %-------------------------------------------------------------------------
                 %Obtain the noise
    
                 sigma_Z = 0;    % noise variance
                 w= get_noise(m,sigma_Z);
                 
                %-------------------------------------------------------------------------
                % Obtain the compressible signal

                %epsilon = 1e-15;
                sigma_1 = 10;   % variance of distribution representing large coefficients
                sigma_0 = 1;    % variance of distribution representing small coefficients 

                [x supp_index partial_supp_index cardinal_par_supp_set] = get_compressible_signal(n,k,sigma_1,sigma_0);
                
                support_set=zeros(1,n);
                support_set(supp_index)=1;
                partial_supp_set=ones(1,n)*((k-cardinal_par_supp_set)/(n-cardinal_par_supp_set));
                partial_supp_set(partial_supp_index)= 1;
                

                %-------------------------------------------------------------------------
                %  Random Gaussian _ encoding and decoding
                %-------------------------------------------------------------------------
                
               
                %-------------------------------------------------------------------------
                %  Obtain the Random Gaussian encoding matrix 

                A = encode_random_gaussian_matrix(m,n);
                A = 1/sqrt(m)*A;   
                

                %-------------------------------------------------------------------------
                %Noisy Measurement 

                y = A*x + w;
                
                save('data_sample_500');