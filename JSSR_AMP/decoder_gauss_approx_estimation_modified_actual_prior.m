function [xrec,q]=decoder_gauss_approx_estimation_modified_actual_prior(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,num_outer_iter,damp_factor_value_xrec,damp_factor_q)
%---------------------------------------------------------
% initializations
%---------------------------------------------------------
%    clc;
%    clear all;
%   %load 'data_sample_500';
% [x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter, num_outer_iter , support_set] = get_input_data();
 num_inner_iter =150;
xrec=zeros(n,num_inner_iter,1); % max likelihood
damp_factor_value_xrec = 0.3;
damp_factor_q = 0.9;
s=k/n;
% epsilon=1e-15;
% sigma_0 =epsilon;
q=repmat(s,num_inner_iter,n);
variance_1 = sigma_1*sigma_1;
variance_0 = sigma_0*sigma_0;
mean_1 = 0;
mean_0 = 0;



   %------------------------------------------------------------------   
   % BP iterations
   %------------------------------------------------------------------ 
    
   
         
   %------------------------------------------------------------------   
   % updates requied for the inner iterations
   %------------------------------------------------------------------ 
   
   
   
   brdcast_message_var_to_constr_weight = repmat(s,1,n);
   
   approx_prior_variance = s*variance_1 + (1-s)*variance_0 + s.*(1-s).*((mean_1 - mean_0).^2);
   brdcast_message_var_to_constr_weighted_mean_sum = repmat(0,1,n);
   brdcast_message_var_to_constr_inverse_variance_sum =approx_prior_variance;
             
       
   for iter_count=1:num_inner_iter	 
       disp('iter_conunt=');
       disp(iter_count);
      
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
                   ret_prod_incoming_messages_variance =  1./ret_prod_incoming_messages_inv_variance;             
                   
                   ret_prod_incoming_message_mean = (repmat(brdcast_message_var_to_constr_weighted_mean_sum,m,1)- ...
                       message_constr_to_var_node_mean.*message_constr_to_var_node_variance_inverse)...
                    .*ret_prod_incoming_messages_variance;
           

           rep_var_0 =  repmat(variance_0,m,n);
           rep_mean_0 = repmat(mean_0,m,n);
           rep_var_1 = repmat(variance_1,m,n);
           rep_mean_1 = repmat(mean_1,m,n);
           rep_weight = repmat(q(iter_count-1,:),m,1);
           
           [weight_0_new mean_final_0 variance_final_0] = weight_determination(rep_var_0,rep_mean_0,ret_prod_incoming_messages_variance,ret_prod_incoming_message_mean,1-rep_weight);
           [weight_1_new mean_final_1 variance_final_1] = weight_determination(rep_var_1,rep_mean_1,ret_prod_incoming_messages_variance,ret_prod_incoming_message_mean,rep_weight);
           
            small_index_weight_0 =  find(weight_0_new<=1e-30);
            weight_0_new(small_index_weight_0) = 0;
            weight_1_new(small_index_weight_0)= 1;
           
            weight_0_final = weight_0_new ./ (weight_0_new + weight_1_new);
            weight_1_final = weight_1_new ./ (weight_0_new + weight_1_new);
 
            
            
            approx_mean = weight_0_final.*mean_final_0 + weight_1_final.*mean_final_1;
            approx_variance = weight_0_final.*variance_final_0 + weight_1_final.*variance_final_1 +...
                weight_0_final.*weight_1_final.*(mean_final_0 - mean_final_1).^2;
            
           damp_factor_ret_mean =0.5;     
           ret_message_var_to_constr_mean =  damp_factor_ret_mean * ret_message_var_to_constr_mean  + (1-damp_factor_ret_mean)*approx_mean;
                  
           ret_message_var_to_constr_variance = approx_variance;
            
           end
   
    
           %------------------------------------------------------------------   
           % Convolve the pdfs and update the messages
           %------------------------------------------------------------------
            
            A_square = A.*A;
            variance_conv_temp = ret_message_var_to_constr_variance.*A_square;
            message_constr_to_var_node_variance = (repmat([sum(variance_conv_temp')]',1,n) - variance_conv_temp)./A_square;
            mean_conv_temp = ret_message_var_to_constr_mean.*A;
            message_constr_to_var_node_mean = (repmat([sum(mean_conv_temp')]',1,n) - mean_conv_temp);
            message_constr_to_var_node_mean = (repmat(y,1,n) - message_constr_to_var_node_mean)./A;
      
	
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
         % Take product of all the incoming messages and send out the
         % weighted sum of means and variances  and  send it out 
         %-----------------------------------------------------------------
         
         
           message_constr_to_var_node_variance_inverse = 1./message_constr_to_var_node_variance;
           brdcast_message_var_to_constr_inverse_variance_sum = sum(message_constr_to_var_node_variance_inverse);
           
           brdcast_message_var_to_constr_weighted_mean_sum = sum([message_constr_to_var_node_mean.*message_constr_to_var_node_variance_inverse]);
           
         %-----------------------------------------------------------------  
         % 1.a.Calculate the belief
         %-----------------------------------------------------------------
           prod_variance = 1./brdcast_message_var_to_constr_inverse_variance_sum; % variance of the product of all the incoming message densities
           prod_mean = brdcast_message_var_to_constr_weighted_mean_sum.*prod_variance; % mean of the product of all the incoming message densities
           
           rep_var_0 =  repmat(variance_0,1,n);
           rep_mean_0 = repmat(mean_0,1,n);
           rep_var_1 = repmat(variance_1,1,n);
           rep_mean_1 = repmat(mean_1,1,n);
           rep_weight = q(iter_count,:);
           
           [weight_0_new mean_final_0 variance_final_0] = weight_determination(rep_var_0,rep_mean_0,prod_variance,prod_mean,1-rep_weight);
           [weight_1_new mean_final_1 variance_final_1] = weight_determination(rep_var_1,rep_mean_1,prod_variance,prod_mean,rep_weight);
                
         
         %------------------------------------------------------------------   
         % 3.Calculate the estimate of the signal
         %------------------------------------------------------------------
       for var_node_cnt = 1:n  
         pdf_max_value_1 = weight_0_new(var_node_cnt) * normpdf(mean_final_0(var_node_cnt),...
                         mean_final_0(var_node_cnt),sqrt(variance_final_0(var_node_cnt))) + ...
                         weight_1_new(var_node_cnt) * normpdf(mean_final_0(var_node_cnt),...
                         mean_final_1(var_node_cnt),sqrt(variance_final_1(var_node_cnt))); 
         pdf_max_value_2 = weight_1_new(var_node_cnt) * normpdf(mean_final_1(var_node_cnt),...
                         mean_final_1(var_node_cnt),sqrt(variance_final_1(var_node_cnt))) + ...
                         weight_0_new(var_node_cnt) * normpdf(mean_final_1(var_node_cnt),...
                         mean_final_0(var_node_cnt),sqrt(variance_final_0(var_node_cnt))); 
                     
         if (pdf_max_value_1 > pdf_max_value_2)
             xrec_temp(var_node_cnt)= mean_final_0(var_node_cnt);
         else
             xrec_temp(var_node_cnt)= mean_final_1(var_node_cnt);
         end        
       end        
        
       
            %damp_factor_value_xrec = min(abs(10./xrec_temp'),1);
            
            
            xrec(:,iter_count) = damp_factor_value_xrec.*xrec(:,max(iter_count-1,1)) + (1-damp_factor_value_xrec).*xrec_temp';
         %   test_xrec(:,iter_count) = xrec(:,outer_iter_count);
           
	           

     
      %------------------------------------------------------------------   
     % Update the value of Q from estimate of x
     %------------------------------------------------------------------
         
         prob_q_large = [normpdf(xrec(:,iter_count),0,sigma_1)]'.*q(iter_count,:);
         prob_q_small = [normpdf(xrec(:,iter_count),0,sigma_0)]'.*[1-q(iter_count,:)];
         q_sum = prob_q_large+prob_q_small;
         
         small_index = find( q_sum < 1e-30);% if the estimate is too large declare that component as large 
         prob_q_large(small_index) = 0.8;
         prob_q_small(small_index) = 0.2;
         temp_q = prob_q_large./(prob_q_large + prob_q_small);
         damp_factor_q=0.9;
%          q_threshold = 0.7;
%          thresholded_q = (sign(init_esti_q - q_threshold)+1)/2 ;
%          temp_q = min( thresholded_q + (k - sum(thresholded_q))/n,1);    
         q(iter_count+1,:)=damp_factor_q*q(iter_count,:) + (1-damp_factor_q)*temp_q;
         
   end %of inner iterations            
     
    disp('end of decoder unknown support')
    


%---------------------------------------------
% checked on 15th dec
% dec 20 - check for the variance equation