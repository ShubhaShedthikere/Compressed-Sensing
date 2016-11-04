function [xrec,q]=decoder_gauss_approx_estimation_modified_approx_prior(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,damp_factor_q)
%---------------------------------------------------------
% initializations
%---------------------------------------------------------
%     clc;
%     clear all;
%    %load 'data_sample_500';
  load 'slopes_threshold_sigma_0_large';
%   [x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter, num_outer_iter , support_set] = get_input_data();
 num_inner_iter =150;
xrec=zeros(n,num_inner_iter,1); % max likelihood
s=k/n;
q=repmat(s,num_inner_iter,n);
% epsilon=1e-15;
% sigma_0 =epsilon;
variance_1 = sigma_1*sigma_1;
variance_0 = sigma_0*sigma_0;
mean_1 = 0;
mean_0 = 0;
A_square = A.*A;

message_constr_to_var_node_mean = zeros(m,n);


%tic
   %------------------------------------------------------------------   
   % BP iterations
   %------------------------------------------------------------------ 
    
   
         
   %------------------------------------------------------------------   
   % updates requied for the inner iterations
   %------------------------------------------------------------------ 
     
   for iter_count=1:num_inner_iter	  
      
    %------------------------------------------------------------------   
    % Processing at the constraint node 
    % ---------------------------------
    % 1.Retrive the correct message from the broadcast message
    %   ( no retrival required for the first iteration as the incoming
    %   messages are the priors)
    % 2.Approximate the gaussian mixtures to single mixture
    % 3.Convolve the incoming pdfs 
    % 4.Update outgoing messages  
    %------------------------------------------------------------------
  
           %------------------------------------------------------------------   
           % Extract the mean and variance and find the approximate
           % gaussian function
           %------------------------------------------------------------------
           if( iter_count == 1)
                    avg_signal_energy = q(iter_count,:)* variance_1 + (1-q(iter_count,:))*variance_0;
                    ret_message_var_to_constr_variance = repmat(avg_signal_energy,m,1);
                    ret_message_var_to_constr_mean = repmat(0,m,n);
                    
           else
                   ret_prod_incoming_messages_inv_variance = repmat(brdcast_message_var_to_constr_inverse_variance_sum,m,1)-...
                                                                 message_constr_to_var_node_variance_inverse; 
                                     
                   ret_prod_incoming_message_mean_temp = (repmat(brdcast_message_var_to_constr_weighted_mean_sum,m,1)- ...
                       message_constr_to_var_node_mean.*message_constr_to_var_node_variance_inverse);
                             
                   ret_prod_incoming_messages_variance = 1./ret_prod_incoming_messages_inv_variance;
                   ret_prod_incoming_message_mean_temp = ret_prod_incoming_message_mean_temp.*ret_prod_incoming_messages_variance;
                    
                    
           
           damp_factor_ret_mean =0.4;   
%            rep_scale_factor = repmat(1,m,n);
%            ret_message_var_to_constr_mean_scaled = rep_scale_factor.* ret_prod_incoming_message_mean_temp;
           ret_message_var_to_constr_mean = damp_factor_ret_mean*ret_message_var_to_constr_mean + (1-damp_factor_ret_mean)*ret_prod_incoming_message_mean_temp;
                  
           ret_message_var_to_constr_variance = ret_prod_incoming_messages_variance;
            
           end
   
    
           %------------------------------------------------------------------   
           % Convolve the pdfs and update the messages
           %------------------------------------------------------------------            
           
            variance_conv_temp = ret_message_var_to_constr_variance.*A_square;
            message_constr_to_var_node_variance = (repmat([sum(variance_conv_temp')]',1,n) - variance_conv_temp)./A_square;
            
            mean_conv_temp = ret_message_var_to_constr_mean.*A;
            message_constr_to_var_node_mean_temp = (repmat([sum(mean_conv_temp')]',1,n) - mean_conv_temp);
            message_constr_to_var_node_mean_temp = (repmat(y,1,n) - message_constr_to_var_node_mean_temp)./A;
            damp_factor_constr_mean = 0;            
            
            message_constr_to_var_node_mean = message_constr_to_var_node_mean_temp*(1-damp_factor_constr_mean)...
                                              + damp_factor_constr_mean*message_constr_to_var_node_mean;
                                          
	
  %------------------------------------------------------------------   
  % processing at variable node:-
  %-----------------------------
  % 1.a. Calculate the belief 
  % 1.b. Damping the belief update  
  % 2.a. Calculate the outgoing message from the variable node
  % 2.b. Damp the outgoing messages
  % 3.   Find the estimate of the signal
  %------------------------------------------------------------------         
         
         %----------------------------------------------------------------- 
         % Take product of all the incoming messages and find the mean and
         % variance of the incoming messages
         %-----------------------------------------------------------------
                  
           message_constr_to_var_node_variance_inverse = 1./message_constr_to_var_node_variance;
           message_constr_to_var_node_variance_inverse_sum = sum(message_constr_to_var_node_variance_inverse);
           message_var_to_constr_weighted_mean_sum = sum([message_constr_to_var_node_mean.*message_constr_to_var_node_variance_inverse]);
           prod_variance = 1./message_constr_to_var_node_variance_inverse_sum; % variance of the product of all the incoming message densities
           prod_mean = message_var_to_constr_weighted_mean_sum.*prod_variance; % mean of the product of all the incoming message densities
           
           
         %-----------------------------------------------------------------  
         % Calculate the Outgoing messages
         %-----------------------------------------------------------------
           approx_prior_mean=0;            
           approx_prior_variance  = q(iter_count,:)* variance_1 + (1-q(iter_count,:))*variance_0;
           brdcast_message_var_to_constr_inverse_variance_sum = message_constr_to_var_node_variance_inverse_sum + 1./approx_prior_variance;
           brdcast_message_var_to_constr_weighted_mean_sum =  message_var_to_constr_weighted_mean_sum;
           
         %------------------------------------------------------------------   
         % 3.Calculate the final mean and variance after including prior
         %------------------------------------------------------------------  
           
           final_variance = approx_prior_variance.*prod_variance./(approx_prior_variance + prod_variance);
           final_mean =   final_variance .*(prod_mean./prod_variance + approx_prior_mean./approx_prior_variance);
           
         %------------------------------------------------------------------   
         % 4.Determine scaling factor ( this depends on prior)
         %------------------------------------------------------------------  
           mean_threshold_lower = 3; % needs to be tuned!

         %--------------------------------------------------------------------
         %find the slope and threshold
          
         q1= q(iter_count,:);
         
         slope_threshold_lin_amp = interp1(q_table,slopes_threshold_lin_amp_table,q1,'linear','extrap'); % slope is function of q
         threshold_lin_amp = prod_mean.*prod_mean.*slope_threshold_lin_amp;
         slope_lin_amp_region =  interp1(q_table,slopes_lin_amplification_table,q1,'linear','extrap');
    
        %----------------------------------------------------------------------
        %find the region of operation
        %----------------------------------------------------------------------    
        scale_factor = (1+sign(abs(prod_mean) - mean_threshold_lower))/2.*...
                      ((1+sign(prod_variance - threshold_lin_amp))/2.*(0.96*exp(-0.4*prod_variance) + 0.04)... 
                     +(1-sign(prod_variance - threshold_lin_amp))/2.*(1+slope_lin_amp_region.*prod_variance)) + ...
                      (1-sign(abs(prod_mean) - mean_threshold_lower))/2.*(0.96*exp(-0.4*prod_variance) + 0.04);
                         
           
         %------------------------------------------------------------------   
         % 4.Find the estimate of the signal
         %-----------------------------------------------------------------
         
            x_est = final_mean .* scale_factor;
            damp_factor_value_xrec = 0.3;
            
            xrec(:,iter_count) = damp_factor_value_xrec.*xrec(:,max(iter_count-1,1)) + (1-damp_factor_value_xrec).*x_est';      

     
      %------------------------------------------------------------------   
     % Update the value of Q from estimate of x
     %----------------------------------------------------------------
         damp_factor_q = 0;
         prob_q_large = [normpdf(xrec(:,iter_count),0,sigma_1)]'.*q(iter_count,:);
         prob_q_small = [normpdf(xrec(:,iter_count),0,sigma_0)]'.*[1-q(iter_count,:)];
         q_sum = prob_q_large+prob_q_small;
         
         small_index = find( q_sum < 1e-30);% if the estimate is too large declare that component as large 
         prob_q_large(small_index) = 0.8;
         prob_q_small(small_index) = 0.2;
         temp_q = prob_q_large./(prob_q_large + prob_q_small);
         
%          q_threshold = 0.7;
%          thresholded_q = (sign(init_esti_q - q_threshold)+1)/2 ;
%          temp_q = min( thresholded_q + (k - sum(thresholded_q))/n,1);    
         q(iter_count+1,:)=damp_factor_q*q(iter_count,:) + (1-damp_factor_q)*temp_q;
         
       %  disp(iter_count);
   end %of inner iterations            
     
%    diff_approx_prior = repmat(x,1,150)-xrec;
%    mse_approx_prior = sqrt(sum(diff_approx_prior.*diff_approx_prior));
%    plot(mse_approx_prior,'r-');
%    grid on
%      figure()
%      plot(xrec')
%      
%     time_taken= toc;
%     disp('Time Taken : ');
%     disp(time_taken)



%---------------------------------------------
% checked on 15th dec
% dec 20 - check for the variance equation
% 28-Jul-10 : changed the structure for speed.
%             Added scaling factor code
%             itself.