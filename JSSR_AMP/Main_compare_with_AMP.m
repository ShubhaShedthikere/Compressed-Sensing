clear all;
clc;

tic
num_data_sets = 100; % number of realizations


step_size_percentage_sparsity = 0.1;
num_steps_perc_sparsity = 1;

%rho = 0.2;
%delta = 0.3

step_size_n =500;
num_steps_n = 1;

%step_size_m = floor(delta*step_size_n);
step_size_m = 50
num_steps_m = 10;


for spar_cnt= 1:num_steps_perc_sparsity
    
for n_cnt = 1:num_steps_n
    n = n_cnt*step_size_n;
    k = floor(spar_cnt*step_size_percentage_sparsity*n);
    
    %k = floor(rho*step_size_m);
    
    for m_cnt=1:num_steps_m
        m=m_cnt*step_size_m ;
        
    
        
        
                for data_set_cnt=1:num_data_sets
             
                 disp('m =');
                 disp(m)
                 disp('  data_set')
                 disp(data_set_cnt);
                 
                 %-------------------------------------------------------------------------
                 %Obtain the noise
    
                 sigma_Z = 0;    % noise variance
                 w= get_noise(m,sigma_Z);
                 
                %-------------------------------------------------------------------------
                % Obtain the compressible signal

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

                A = get_encoding_matrix_bernoulli(m,n);
                A = 1/sqrt(m)*A;   
                

                %-------------------------------------------------------------------------
                %Noisy Measurement 

                y = A*x + w;
                
                %-------------------------------------------------------------------------
                %Iterative Thresholding Algorithm
                 
                  damp_factor_q = 0.9;
                  num_inner_iter=150;
                  [x_rec_bp_outer_iter_modified_approx_prior, q_bp_outer_iter_modified_approx_prior]  = decoder_gauss_approx_estimation_modified_approx_prior(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,damp_factor_q);
                  error_rec_bp_outer_iter_modified_approx_prior(data_set_cnt,m_cnt) = (norm(x-x_rec_bp_outer_iter_modified_approx_prior(:,num_inner_iter),2));
                  
                  diff_approx_prior = repmat(x,1,num_inner_iter)-x_rec_bp_outer_iter_modified_approx_prior;
                   mse_approx_prior(data_set_cnt,:) = sqrt(sum(diff_approx_prior.*diff_approx_prior))/norm(x,2);
                    
%                    damp_factor_q = 0.9;
%                   num_inner_iter=150;
%                   [x_rec_bp_outer_iter_modified_actual_prior, q_bp_outer_iter_modified_actual_prior]  = decoder_gauss_approx_estimation_modified_actual_prior(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,damp_factor_q);
%                   error_rec_bp_outer_iter_modified_actual_prior(data_set_cnt,m_cnt) = (norm(x-x_rec_bp_outer_iter_modified_actual_prior(:,num_inner_iter),2));
%                   
%                   diff_actual_prior = repmat(x,1,num_inner_iter)-x_rec_bp_outer_iter_modified_actual_prior;
%                    mse_actual_prior(data_set_cnt,:) = sqrt(sum(diff_actual_prior.*diff_actual_prior))/norm(x,2)
%                   
                   
                  damp_factor_q = 0.9;
                  [x_rec_AMP]  = decoder_AMP(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,damp_factor_q);
                  error_rec_AMP(data_set_cnt,m_cnt) = (norm(x-x_rec_AMP(:,num_inner_iter),2));
                  
                  diff_AMP = repmat(x,1,num_inner_iter)-x_rec_AMP;
                    mse_AMP(data_set_cnt,:) = sqrt(sum(diff_AMP.*diff_AMP))/norm(x,2) ;
                  
                
                end        

            
    end% end of the loop for range of m
    
    m = step_size_m*(1:num_steps_m);
    
   % plot(m,avg_error_known_supp,'r-*',m,avg_error_partial_supp(:,1),'b-x',m,avg_error_partial_supp(:,2),'k-x',m,avg_error_partial_supp(:,4),'g-o');
  end% end of loop for range of n

 end % end of loop for range of sparsity


figure()
plot(mse_AMP,'r');
grid on
hold on
plot(mse_approx_prior');



disp('done! ');

toc







