function [xrec]=decoder_AMP(x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter,num_outer_iter,damp_factor_q)
%---------------------------------------------------------
% initializations
%---------------------------------------------------------
%     clc;
%     clear all;
% %    %load 'data_sample_500';
% [x,y,n,k,m,A,sigma_1,sigma_0,sigma_Z,num_inner_iter, num_outer_iter , support_set] = get_input_data();
 num_inner_iter = 150;
xrec=zeros(n,num_inner_iter,1); % max likelihood
s=k/n;
q=repmat(s,num_inner_iter,n);
% epsilon=1e-15;
% sigma_0 =epsilon;
variance_1 = sigma_1*sigma_1;
variance_0 = sigma_0*sigma_0;
mean_1 = 0;
mean_0 = 0;
delta =m/n;



   %------------------------------------------------------------------   
   % BP iterations
   %------------------------------------------------------------------ 
    
   approx_prior_init_variance = s*variance_1 + (1-s)*variance_0 + s.*(1-s).*((mean_1 - mean_0).^2);
   threshold = approx_prior_init_variance *0.25;
   message_var_to_constr_mean = zeros(n,1);     
   message_constr_to_var_node_mean = zeros(m,1);
       
   for iter_count=1:num_inner_iter	  
      
    %------------------------------------------------------------------   
    % Processing at the constraint node 
    % ---------------------------------
    % 1.Convolve the incoming pdfs 
    % 2.Update outgoing messages  
    %------------------------------------------------------------------
  
            
            constr_to_var_node_mean_temp = A*message_var_to_constr_mean;            
            
            constr_to_var_node_mean_temp = y - constr_to_var_node_mean_temp + ...
                                           1/delta*message_constr_to_var_node_mean*sum(sign(abs(message_var_to_constr_mean)))/n;
            
            damp_factor_constr_mean = 0;
            
            message_constr_to_var_node_mean = constr_to_var_node_mean_temp*(1-damp_factor_constr_mean)...
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
         % Take product of all the incoming messages and find the mean 
         %-----------------------------------------------------------------
                  
           prod_variance = threshold/delta;
           
           var_to_constr_mean_temp =  A'*message_constr_to_var_node_mean + message_var_to_constr_mean;
           
           message_var_to_constr_mean = soft_threshold(var_to_constr_mean_temp,prod_variance);
           
           xrec(:,iter_count) = message_var_to_constr_mean;
           
                   
         %-----------------------------------------------------------------  
         % Update threshold
         %-----------------------------------------------------------------
         
           threshold = prod_variance * sum(sign(abs(message_var_to_constr_mean)))/n;
           
     
   end %of inner iterations      
   
%     diff_approx_prior = repmat(x,1,num_inner_iter)-xrec;
%     mse_approx_prior = sqrt(sum(diff_approx_prior.*diff_approx_prior))/norm(x,2);
%     plot(mse_approx_prior,'r');
    
     
    disp('end of decoder unknown support')



%---------------------------------------------
% checked on 15th dec
% dec 20 - check for the variance equation