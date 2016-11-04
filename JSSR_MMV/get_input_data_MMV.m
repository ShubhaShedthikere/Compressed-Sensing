function [x,y,n,k,m,l,A,sigma_1,sigma_s,sigma_Z,num_inner_iter,support_set,w] = get_input_data_MMV()

% rho = 0.3;
% delta = 0.3

n=500;
m=150;
k=50;
l=10;
s=k/n;
num_inner_iter=75;
SNR = 10;

                %-------------------------------------------------------------------------
                % Obtain the compressible signal

                sigma_1 = 10;   % variance of distribution representing large coefficients
                sigma_0 = 0;    % variance of distribution representing small coefficients 

                [x support_set] = get_compressible_signal_MMV(n,k,sigma_1,sigma_0,l);
                       
                %-------------------------------------------------------------------------
                %  Obtain the Random Gaussian encoding matrix 

                A = encode_random_gaussian_matrix(m,n);
                %A_gauss = 1/sqrt(m)*A;   
                

                %-------------------------------------------------------------------------
                %Noiseless measurement

                y_temp = A*x;                
              
                 %-------------------------------------------------------------------------
                 %Obtain the noise whose variance is adjusted so that each
                 %of the L measurements has the same SNR
                 
                 
                 sigma_Z = sqrt((sum(y_temp.*y_temp))/m * 10^(-SNR/10));
                 %sigma_Z=zeros(1,l);
                 
                 w= get_noise(m,sigma_Z,l);
                 
                 %---------------------------------------------------------
                 % noisy measurement
                
                  y= y_temp + w;         
                  
                  sigma_s=1;
               
              
                
                save('data_sample_500_noisy','x','y','n','k','m','l','A','sigma_1','sigma_s','sigma_Z','num_inner_iter','support_set');
                
                %y=y_temp;
                %sigma_Z = 0;
                %estimated_sparsity = ((sum(y.*y)/(m*n)) - sigma_0*sigma_0)/(sigma_1*sigma_1-sigma_0*sigma_0);
                
                %save('data_sample_noisy_l2_test_10db','x','y','n','k','m','l','A','sigma_1','sigma_0','sigma_Z','num_inner_iter','support_set');
                