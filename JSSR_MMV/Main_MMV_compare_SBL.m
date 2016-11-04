clear all;
clc;

tic
num_data_sets = 20; % number of realizations


%step_size_percentage_sparsity = 0.1;
num_steps_perc_sparsity = 1;

n = 500;
%m = 150;
k = 50;
s=k/n;

M_table=[ 100 150 200 250 300 350];

M_num_of_steps = 6;
% L_step_size = 5;
% L_range = 50;
% L_init = 1 ;
% L_num_of_steps = (L_range - L_init)/L_step_size ;
 


L_table = [ 6 7 8 9 10 15 20 ];
L_num_of_steps = 7 ;

data_count=zeros(L_num_of_steps,M_num_of_steps);
      for L_cnt = 1:L_num_of_steps
          
          l = L_table(L_cnt);
    
          for m_cnt = 1:M_num_of_steps;
              
              m = M_table(m_cnt);
              
                 for data_set_cnt=1:num_data_sets
             
                 disp('l =');
                 disp(l)
                 disp('  data_set')
                 disp(data_set_cnt);
                 
                               
                %-------------------------------------------------------------------------
                % Obtain the compressible signal

                sigma_1 = 10;   % variance of distribution representing large coefficients
                sigma_0 = 0;    % variance of distribution representing small coefficients 

                [x support_set] = get_compressible_signal_MMV(n,k,sigma_1,sigma_0,l);
                       
                %-------------------------------------------------------------------------
                %  Obtain the Random Gaussian encoding matrix 

                A = encode_random_gaussian_matrix(m,n);
                %A_gauss = 1/sqrt(m)*A_gauss;   
                

                %-------------------------------------------------------------------------
                %Noiseless measurement

                y_temp = A*x;                
              
                 %-------------------------------------------------------------------------
                 %Obtain the noise whose variance is adjusted so that each
                 %of the L measurements has the same SNR
                 
                 %sigma_Z = zeros(1,l);
                  SNR=10;
                  sigma_Z = sqrt((sum(y_temp.*y_temp))/m * 10^(-SNR/10));
%              
                 w= get_noise(m,sigma_Z,l);
                 
                 %---------------------------------------------------------
                 % noisy measurement
                
                  y= y_temp + w;
               
                  
                 %---------------------------------------------------------
                 %joint sparse decoding known q
                  
                   num_inner_iter=100;
                   sigma_s=1;
                   [x_rec_joint_sparse, q_rec_joint_sparse]  = approx_prior_based_joint_approx_sparse_recovery(x,y,n,k,m,l,A,sigma_1,sigma_s,sigma_Z,num_inner_iter);
                   
                   %find the MSE
                   diff = x-x_rec_joint_sparse;
                   mse_joint_sparse_fro(data_set_cnt,m_cnt,L_cnt) = norm(diff,'fro');
                   mse_fro_temp = mse_joint_sparse_fro(data_set_cnt,m_cnt,L_cnt)^2/norm(x,'fro')^2;
                   mse_temp =  sqrt(mse_joint_sparse_fro(data_set_cnt,m_cnt,L_cnt)*mse_joint_sparse_fro(data_set_cnt,m_cnt,L_cnt)/l);
                   
                   if(mse_temp < 200)
                   mse_joint_sparse_vec_wise(data_set_cnt,m_cnt,L_cnt) = mse_temp;
                   mse_joint_sparse_fro_store(data_set_cnt,m_cnt,L_cnt) = mse_fro_temp;
                   data_count(L_cnt,m_cnt)= data_count(L_cnt,m_cnt)+1;
                   avg_mse_vec_wise(L_cnt,m_cnt) = sum(mse_joint_sparse_vec_wise(:,m_cnt,L_cnt))./data_count(L_cnt,m_cnt)
                   avg_mse_fro(L_cnt,m_cnt) = sum(mse_joint_sparse_fro_store(:,m_cnt,L_cnt))./data_count(L_cnt,m_cnt)                   
                   else
                   mse_joint_sparse_vec_wise(data_set_cnt,m_cnt,L_cnt) = 0;
                   mse_joint_sparse_fro_store(data_set_cnt,m_cnt,L_cnt) = 0;    
                   end
                   
                   
                   %find the error in exact support recovery
                   
                   ber(data_set_cnt,m_cnt,L_cnt) = sum(abs(support_set - q_rec_joint_sparse));
                   
                    save('MMV_noisy_JSDMP_3','mse_joint_sparse_vec_wise','ber','mse_joint_sparse_fro_store','data_count','avg_mse_vec_wise','avg_mse_fro');
                 %---------------------------------------------------------
                 %SBL
                  
                  [x_rec_SBL, q_SBL]  = M_Sparse_Bayesian_Learning(y,A,m,n,sigma_Z,l);
                   
                   %find the MSE
                   diff = x-x_rec_SBL;
                   mse_joint_sparse_fro_SBL(data_set_cnt,m_cnt,L_cnt) = norm(diff,'fro');
                   mse_joint_sparse_fro_store_SBL(data_set_cnt,m_cnt,L_cnt)=mse_joint_sparse_fro_SBL(data_set_cnt,m_cnt,L_cnt)^2/norm(x,'fro')^2;
                   mse_joint_sparse_vec_wise_SBL(data_set_cnt,m_cnt,L_cnt) = sqrt(mse_joint_sparse_fro_SBL(data_set_cnt,m_cnt,L_cnt)*mse_joint_sparse_fro_SBL(data_set_cnt,m_cnt,L_cnt)/l);
                    avg_mse_vec_wise_SBL(L_cnt,m_cnt) = sum(mse_joint_sparse_vec_wise_SBL(:,m_cnt,L_cnt))./data_count(L_cnt,m_cnt)
                   avg_mse_fro_SBL(L_cnt,m_cnt) = sum(mse_joint_sparse_fro_store_SBL(:,m_cnt,L_cnt))./data_count(L_cnt,m_cnt)  
                   
                   %find the error in exact support recovery
                   
                   ber_SBL(data_set_cnt,m_cnt,L_cnt) = sum(abs(support_set - q_SBL')); 
                     save('MMV_noisy_SBL_3','mse_joint_sparse_vec_wise_SBL','ber_SBL','mse_joint_sparse_fro_store_SBL','data_count','avg_mse_vec_wise_SBL','avg_mse_fro_SBL');
%                    
                end % end of datasets loop       
          end % end of m
            
        end% end of the loop for range of L
    



disp('done! ');

toc








