function [xrec_final,supp_rec]=approx_prior_based_joint_approx_sparse_recovery(x,y,n,k,m,l,A,sigma_1,sigma_0,sigma_Z,num_inner_iter)
%---------------------------------------------------------
% initializations
%---------------------------------------------------------
  %clc;
  %clear all;
 %load 'data_sample_500_noisy';
  load 'slopes_threshold_sigma_0_large';
  sigma_0=1;
%[x,y,n,k,m,l,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,support_set] = get_input_data_MMV();
 num_inner_iter=100;

xrec=zeros(n,num_inner_iter,l); % max likelihood
xrec_final = zeros(n,l);
%estimated_sparsity = ((sum(sum(y.*y))/(m*n*l)) - sigma_0*sigma_0)/(sigma_1*sigma_1-sigma_0*sigma_0);
s=k/n;
q=repmat(s,[l,n,num_inner_iter]);
q_given_y(1,:)=repmat(s,1,n);
epsilon=1e-15;

sigma_0 =max(sigma_0,epsilon);
sigma_Z = max(sigma_Z,epsilon);

variance_1 = sigma_1*sigma_1;
variance_0 = sigma_0*sigma_0;
variance_Z = sigma_Z.*sigma_Z;
mean_1 = 0;
mean_0 = 0;
A_square = A.*A;

message_constr_to_var_node_mean = zeros(m,n,l);
message_sum_LLR = zeros(num_inner_iter,n);
message_noise_to_constr_mean = zeros(m,l);

   %------------------------------------------------------------------   
   % BP iterations
   %------------------------------------------------------------------ 
   
  LLR_q_prior = log((s)/(1-s));  % LLR = Log ( p( q=1) / p(q=0))
  
       
  for outer_iter_count =1:1:2
     
message_constr_to_var_node_mean = zeros(m,n,l);
message_sum_LLR = zeros(num_inner_iter,n);
message_noise_to_constr_mean = zeros(m,l);
xrec=zeros(n,num_inner_iter,l); % max likelihood
xrec_final = zeros(n,l);

if(outer_iter_count == 2)
%     sigma_0 = 1e-15;
num_inner_iter =20;
end
q=repmat(s,[l,n,num_inner_iter]);
      
   for iter_count=1:num_inner_iter	  
       
       for mmv_cnt=1:l
      
    %------------------------------------------------------------------------------------------------------------------------------------   
    % Processing at the constraint node 
    %------------------------------------------------------------------------------------------------------------------------------------   
    % 1.Retrive the correct message from the broadcast message
    %   ( no retrival required for the first iteration as the incoming
    %   messages are the priors)
    % 2.Approximate the gaussian mixtures to single mixture
    % 3.Convolve the incoming pdfs 
    % 4.Update outgoing messages  
    %------------------------------------------------------------------
  
           %------------------------------------------------------------------   
           % Extract the correct q,mean and variance from the broadcast
           % message
           %------------------------------------------------------------------
           
           if( iter_count == 1)
                    if(outer_iter_count ==1)
                    q(mmv_cnt,:,iter_count)=repmat(s,1,n);
                    else
                    q(mmv_cnt,:,iter_count)= supp_rec;
                    end
                    avg_signal_energy = q(mmv_cnt,:,iter_count)* variance_1 + (1-q(mmv_cnt,:,iter_count))*variance_0;
                    ret_message_var_to_constr_variance = repmat(avg_signal_energy,m,1);
                    ret_message_var_to_constr_mean = repmat(0,m,n);
                    message_noise_to_constr_variance(:,mmv_cnt) = repmat(variance_Z(mmv_cnt),m,1);
                    message_noise_to_constr_mean(:,mmv_cnt)=zeros(m,1);
                    
           else
                 if(outer_iter_count ==1)  
                   if(l>1)
                   ret_LLR = message_sum_LLR(iter_count-1,:) - message_LLR_q(mmv_cnt,:,iter_count-1);
                   else
                   ret_LLR = message_sum_LLR(iter_count-1,:);   
                   end
                   ret_LR = exp(sign(ret_LLR).*min(abs(ret_LLR),500));
                   temp_q = ret_LR./(ret_LR + 1 );
                   damp_factor_q = 0.9;
                   q(mmv_cnt,:,iter_count)=damp_factor_q*q(mmv_cnt,:,iter_count-1) + (1-damp_factor_q)*temp_q;
                   %LLR_q_prior(iter_count,:) = log(q(mmv_cnt,:,iter_count)./ (1-q(mmv_cnt,:,iter_count)));
                 else
                    q(mmv_cnt,:,iter_count)= supp_rec;
                 end
             
                   constr_to_var_node_variance_inverse_temp = 1./message_constr_to_var_node_variance(:,:,mmv_cnt);
                   ret_prod_incoming_messages_inv_variance = repmat(brdcast_message_var_to_constr_inverse_variance_sum(mmv_cnt,:),m,1)-...
                                                                 constr_to_var_node_variance_inverse_temp; 
                             
                   
                   ret_prod_incoming_message_mean_temp = (repmat(brdcast_message_var_to_constr_weighted_mean_sum(mmv_cnt,:),m,1)- ...
                       message_constr_to_var_node_mean(:,:,mmv_cnt).*constr_to_var_node_variance_inverse_temp);

                   ret_prod_incoming_messages_variance = 1./ret_prod_incoming_messages_inv_variance;
                   
                   ret_prod_incoming_message_mean_temp = ret_prod_incoming_message_mean_temp.*ret_prod_incoming_messages_variance;
                    
                    
           if(outer_iter_count ==1)
           damp_factor_ret_mean =0.5;  %0 for m=50 n==500 k=50
           else
           damp_factor_ret_mean = 0;
           end
%            rep_scale_factor = repmat(1,m,n);
%            ret_message_var_to_constr_mean_scaled = rep_scale_factor.* ret_prod_incoming_message_mean_temp;
           ret_message_var_to_constr_mean = damp_factor_ret_mean*ret_message_var_to_constr_mean + (1-damp_factor_ret_mean)*ret_prod_incoming_message_mean_temp;
                  
           ret_message_var_to_constr_variance = ret_prod_incoming_messages_variance;
            
           end
   
    
           %------------------------------------------------------------------   
           % Convolve the pdfs and update the messages
           %------------------------------------------------------------------
            
            variance_conv_temp = ret_message_var_to_constr_variance.*A_square;
            message_constr_to_noise_variance(:,mmv_cnt) = [sum(variance_conv_temp')]';
            message_constr_to_var_node_variance(:,:,mmv_cnt) = (repmat(message_constr_to_noise_variance(:,mmv_cnt)+ message_noise_to_constr_variance(:,mmv_cnt),1,n) - variance_conv_temp)./A_square;
            
            mean_conv_temp = ret_message_var_to_constr_mean.*A;
            message_constr_to_noise_mean(:,mmv_cnt)= y(:,mmv_cnt)-[sum(mean_conv_temp')]';
            
            constr_to_var_node_mean_temp = (repmat(message_constr_to_noise_mean(:,mmv_cnt) - message_noise_to_constr_mean(:,mmv_cnt),1,n) + mean_conv_temp)./A;
           
            damp_factor_constr_mean = 0;
            
            message_constr_to_var_node_mean(:,:,mmv_cnt) = constr_to_var_node_mean_temp*(1-damp_factor_constr_mean)...
                                              + damp_factor_constr_mean*message_constr_to_var_node_mean(:,:,mmv_cnt);
                                          
	
  %------------------------------------------------------------------------------------------------------------------------------------   
  % Processing at the variable node 
  %------------------------------------------------------------------------------------------------------------------------------------   
  % 1.Extract the q for the particular node from the LLR broadcast information received from the support node
  % 2.Find the mean and variance of product of the messages coming from al
  % the constraint nodes 
  % 3. Use the new q to find the variance of the approximate prior density
  % 4. Add this to the broadcast variance message
  % 5. estimate the value of x using the new q , that is, in finding the
  % scale factor
  % 6. from the estimate of x, find the estimate of q and then LLR
  % 7. Send the broadcast message of mean and variance to constraint node
  % and LLR to the support node
  %------------------------------------------------------------------

         
         %----------------------------------------------------------------- 
         % Calculate the new variance of the prior distribution of x based
         % on the new q value
         %-----------------------------------------------------------------
  
           approx_prior_mean=0;            
           approx_prior_variance  = q(mmv_cnt,:,iter_count)* variance_1 + (1-q(mmv_cnt,:,iter_count))*variance_0;
         
         %----------------------------------------------------------------- 
         % Take product of all the incoming messages and find the mean and
         % variance of the incoming messages
         %-----------------------------------------------------------------
                  
           constr_to_var_node_variance_inverse = 1./message_constr_to_var_node_variance(:,:,mmv_cnt);
           constr_to_var_node_variance_inverse_sum = sum(constr_to_var_node_variance_inverse);
           var_to_constr_weighted_mean_sum = sum([message_constr_to_var_node_mean(:,:,mmv_cnt).*constr_to_var_node_variance_inverse]);
           
           prod_variance = 1./constr_to_var_node_variance_inverse_sum; % variance of the product of all the incoming message densities
           prod_mean = var_to_constr_weighted_mean_sum.*prod_variance; % mean of the product of all the incoming message densities
           
           
         %-----------------------------------------------------------------  
         % Calculate the Outgoing broadcast messages
         %-----------------------------------------------------------------           
           
           brdcast_message_var_to_constr_inverse_variance_sum(mmv_cnt,:) = constr_to_var_node_variance_inverse_sum + 1./approx_prior_variance ;
           brdcast_message_var_to_constr_weighted_mean_sum(mmv_cnt,:) =  var_to_constr_weighted_mean_sum; % anyway the prior mean is 0
           
           
         %------------------------------------------------------------------   
         % Calculate the estimate of the signal
         %------------------------------------------------------------------  
           

           final_variance = approx_prior_variance.*prod_variance./(approx_prior_variance + prod_variance);
           final_mean =   final_variance .*(prod_mean./prod_variance + approx_prior_mean./approx_prior_variance);                 
             
          if(outer_iter_count ==1)
              
         %------------------------------------------------------------------   
         % 4.Determine scaling factor ( this depends on prior)
         %------------------------------------------------------------------  
           mean_threshold_lower = 1; % needs to be tuned! 3 for sigma_0 = 1 0 for sigma_0 = 1e-15 1 for sigma_1=0.3162

         %-----------------------------------------------------------------
         %find the slope and threshold
          
         q1=q(mmv_cnt,:,iter_count);
         
         slope_threshold_lin_amp = interp1(q_table,slopes_threshold_lin_amp_table,q1,'linear','extrap'); % slope is function of q
         threshold_lin_amp = prod_mean.*prod_mean.*slope_threshold_lin_amp;
         slope_lin_amp_region =  interp1(q_table,slopes_lin_amplification_table,q1,'linear','extrap');
    
        %----------------------------------------------------------------------
        %find the region of operation
        %----------------------------------------------------------------------    
        scale_factor = (1+sign(abs(prod_mean) - mean_threshold_lower))/2.*...
                      ((1+sign(prod_variance - threshold_lin_amp))/2.*( 0.96*exp(-0.4*prod_variance) + 0.04)... 
                     +(1-sign(prod_variance - threshold_lin_amp))/2.*(1+slope_lin_amp_region.*prod_variance)) + ...
                      (1-sign(abs(prod_mean) - mean_threshold_lower))/2.*( 0.96*exp(-0.4*prod_variance) + 0.04);
          else
              scale_factor=1;
          end
           
         %------------------------------------------------------------------   
         % 4.Find the estimate of the signal
         %-----------------------------------------------------------------
         
           x_est = final_mean .* scale_factor;
           damp_factor_value_xrec = 0.3;
            
           xrec(:,iter_count,mmv_cnt) = damp_factor_value_xrec.*xrec(:,max(iter_count-1,1),mmv_cnt) + (1-damp_factor_value_xrec).*x_est';    
           xrec_final(:,mmv_cnt) = damp_factor_value_xrec.*xrec_final(:,mmv_cnt) + (1-damp_factor_value_xrec).*x_est'; 

           
       if(outer_iter_count ==1)
     
        %------------------------------------------------------------------   
        % Update the LLR to be passed to the support nodes
        %----------------------------------------------------------------
                        
                         
         prob_q_large =  exp(-0.5 * ((xrec(:,iter_count,mmv_cnt) - 0)./sigma_1).^2) ./ (sqrt(2*pi) .* sigma_1);
         prob_q_small =  exp(-0.5 * ((xrec(:,iter_count,mmv_cnt) - 0)./sigma_0).^2) ./ (sqrt(2*pi) .* sigma_0);
         q_sum = prob_q_large+prob_q_small;
         
         small_index = find( prob_q_small < 1e-300);% if the estimate is too large declare that component as large 
        %prob_q_large(small_index) = 0.8;
         prob_q_small(small_index) = 1e-300;
         
         message_LLR_q(mmv_cnt,:,iter_count) = log(max(prob_q_large./prob_q_small,1e-300));
         
         
       %------------------------------------------------------------------------------------------------------------------------------------   
       % Processing at the noise node 
       %------------------------------------------------------------------------------------------------------------------------------------   
         
%          damp_noise_mean=0;
%          noise_prior = variance_Z(mmv_cnt);
%          message_noise_to_constr_variance(:,mmv_cnt) = (message_constr_to_noise_variance(:,mmv_cnt)*noise_prior)./(message_constr_to_noise_variance(:,mmv_cnt)+noise_prior);
%          message_noise_to_constr_mean_temp = (message_constr_to_noise_mean(:,mmv_cnt)./message_constr_to_noise_variance(:,mmv_cnt)).*message_noise_to_constr_variance(:,mmv_cnt);
%          message_noise_to_constr_mean(:,mmv_cnt) = message_noise_to_constr_mean_temp*(1-damp_noise_mean) + message_noise_to_constr_mean(:,mmv_cnt)*damp_noise_mean;
%   
           end
         %disp(iter_count);
         end % for all signal vectors  
         
        %------------------------------------------------------------------------------------------------------------------------------------   
        % Processing at the support node 
        %------------------------------------------------------------------------------------------------------------------------------------   
        
        if(outer_iter_count==1)
         damp_factor_LLR_q = 0;
         if(l>1)
         message_sum_LLR_temp = sum(message_LLR_q(:,:,iter_count)) + LLR_q_prior;
         else
         message_sum_LLR_temp =    message_LLR_q(:,:,iter_count) + LLR_q_prior;
         end
         
         message_sum_LLR(iter_count,:) = message_sum_LLR(max(iter_count-1,1),:)*damp_factor_LLR_q + (1-damp_factor_LLR_q)*message_sum_LLR_temp;
         
         LR_q_temp = exp(sign(message_sum_LLR(iter_count,:)).*min(abs(message_sum_LLR(iter_count,:)),500));
         
         damp_factor_q_given_y = 0.9;
         q_given_y(iter_count+1,:) = (LR_q_temp./(LR_q_temp + 1 ))*(1-damp_factor_q_given_y) + damp_factor_q_given_y*q_given_y(iter_count,:);
         support_recovered(iter_count,:) = q_given_y(iter_count+1,:);
         
         %LLR_q_prior(iter_count+1,:) = log( q_given_y(iter_count+1,:)./(1-q_given_y(iter_count+1,:)));
          
         diff(iter_count) = norm(y-A*xrec_final,2);
        end
         
   end
   
       
     supp_rec = (sign(support_recovered(num_inner_iter,:)-0.5)+1)/2;
    
     
     
   %diff_approx_prior = x-xrec_final;
  % mse_approx_prior(outer_iter_count,:) = sqrt(sum(diff_approx_prior.*diff_approx_prior))
   
   end % end of outermost for loop
     %err_supp = sum(abs(supp_rec-support_set));
     
     %sprintf('error in support =%d \n MSE = ',err_supp)
    
  % xrec_final = xrec_final.*repmat(supp_rec',1,l);
    disp('end of decoder unknown support')



%---------------------------------------------
% checked on 15th dec
% dec 20 - check for the variance equation
% 28-Jul-10 : changed the structure for speed.
%             Added scaling factor code
%             itself.