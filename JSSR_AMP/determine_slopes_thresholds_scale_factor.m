%function [threshold_lin_amp_table threshold_zero_force_table slopes_threshold_lin_amp_table slopes_threshold_zero_force_table slopes_lin_amplification_table q_table] = determine_slopes_thresholds_scale_factor(sigma_0,sigma_1)

clc
clear all;

load 'ratio_maxactual_max_approx_large_full';
q= q_table;
ratio_size = size(ratio_mean);
num_of_steps_variance = range_upper_variance/variance_step_size;
num_of_steps_mean_prod = range_upper_mean_prod/mean_prod_step_size;
large_number = 10000;
sigma_0=1e-15;
sigma_1 = 10;


for i=1:ratio_size(3)
    
ratio = ratio_mean(:,:,i);
%------------------------------------------------------------------------
%To find the threshold of linear amplification region 
%------------------------------------------------------------------------

threshold_lin_amp_min_diff =0.4;

 
 diff = ratio([1:num_of_steps_variance-1],:)  - ratio([2:num_of_steps_variance],:);
 [max_value max_ind_diff] = max(diff);
 

 
 for temp_cnt=1:1:num_of_steps_mean_prod
     
     if(max_value(temp_cnt) >= threshold_lin_amp_min_diff)
         
          threshold_lin_amp(i,temp_cnt) = max_ind_diff(temp_cnt)*variance_step_size;
          value_at_threshold(i,temp_cnt) =  ratio(max_ind_diff(temp_cnt),temp_cnt);
     else
          
          threshold_lin_amp(i,temp_cnt) = 0;% max(max_ind(temp_cnt)*variance_step_size;
     end
     
 end
 
%------------------------------------------------------------------------
% find the slope by taking the average slope for linear amp region
%------------------------------------------------------------------------
 
 diff_threshold_lin_amp  = (threshold_lin_amp(i,[2:num_of_steps_mean_prod]) - threshold_lin_amp(i,[1:num_of_steps_mean_prod-1]));
 ind_temp = find(diff_threshold_lin_amp >0);
 size_ind_temp = size(ind_temp);
 
 min_mean_value = ind_temp(1)*mean_prod_step_size;
 max_mean_value = (ind_temp(size_ind_temp(2))+1)*mean_prod_step_size; 
 
 slopes_threshold_lin_amp(i) = sum(diff_threshold_lin_amp(ind_temp))/(max_mean_value - min_mean_value);
 
%------------------------------------------------------------------------
% To find the slopes of the amplification factors
%------------------------------------------------------------------------
 lower_end_ind = find(value_at_threshold(i,:)>=1)
 
 lower_end_threshold = threshold_lin_amp(i,lower_end_ind(1));
 lower_end_value = value_at_threshold(i,lower_end_ind(1));
 
 [upper_end_threshold upper_threshold_index]= max(threshold_lin_amp(i,:));
 upper_end_value = value_at_threshold(i,upper_threshold_index);
 
 slopes_lin_amplification(i) = (upper_end_value - lower_end_value)/(upper_end_threshold - lower_end_threshold);
 
%------------------------------------------------------------------------
% To find the threshold zero forcing
%------------------------------------------------------------------------
[min_values min_index] =  min(ratio);

threshold_zero_force(i,:) = min_index*variance_step_size;

epsilon =1e-15;
for temp_cnt=1:1:num_of_steps_mean_prod
     if(min_values(temp_cnt) <= epsilon)
          threshold_zero_force(i,temp_cnt) = min_index(temp_cnt)*variance_step_size;
     else
          threshold_zero_force(i,temp_cnt) = large_number;
     end
     
end

%------------------------------------------------------------------------
% find the slope by taking the average slope for zero force
%------------------------------------------------------------------------
    ind =find(threshold_zero_force(i,:)~=large_number);
    size_ind = size(ind);
    diff_threshold_zero_force = (threshold_zero_force(i,[2:size_ind(2)]) - threshold_zero_force(i,[1:size_ind(2)-1]));
    slopes_threshold_zero_force(i) = sum(diff_threshold_zero_force)./(size_ind(2)-1);

end 

 %------------------------------------------------------------------------
 % Plots
 %------------------------------------------------------------------------
 
 tlt_sigma_0 = num2str(sigma_0);
 tlt_sigma_1 = num2str(sigma_1);
 
 figure(1)
 plot(q,slopes_threshold_zero_force,'-*')
 hold on
 
 figure(2)
 plot(q,slopes_threshold_lin_amp,'-*')
 hold on
 
 figure(3)
 plot(q,slopes_lin_amplification,'-*')
 hold on
  
 figure(1);
 grid on
 xlabel('q')
 ylabel('Slope of the threshold for the zero force region')
 tlt = strcat('Slope of the threshold for the zero force region Vs q for ' ,' Sigma_0 =',tlt_sigma_0,' and Sigma_1=',tlt_sigma_1);
 title(tlt);
 %legend('sigma_0 = 1','sigma_0 = 0.1','sigma_0 = 0.01','sigma_0 = 0.001');
 hgsave('slopes_threshold_zero_force_vs_q_sigma_0_large_init');
  
 figure(2);
 grid on
 xlabel('q')
 ylabel('Slope of the threshold for the linear amplification region')
 tlt = strcat('Slope of the threshold for the linear amplification region Vs q for ' ,' Sigma_0 =',tlt_sigma_0,' and Sigma_1=',tlt_sigma_1);
 title(tlt);
 %legend('sigma_0 = 1','sigma_0 = 0.1','sigma_0 = 0.01','sigma_0 = 0.001');
 hgsave('slopes_threshold_lin_amp_vs_q_sigma_0_large_init');

 figure(3);
 grid on
 xlabel('q')
 ylabel('Slope of the linear amplification factors')
 tlt = strcat('Slope of the linear amplification factors Vs q for ' ,' Sigma_0 =',tlt_sigma_0,' and Sigma_1=',tlt_sigma_1);
 title(tlt);
 %legend('sigma_0 = 1','sigma_0 = 0.1','sigma_0 = 0.01','sigma_0 = 0.001');
 hgsave('slopes_lin_amp_factor_vs_q_sigma_0_large_init');

q_table=q;
threshold_lin_amp_table = threshold_lin_amp;
threshold_zero_force_table = threshold_zero_force;
slopes_threshold_lin_amp_table = slopes_threshold_lin_amp;
slopes_threshold_zero_force_table = slopes_threshold_zero_force;
slopes_lin_amplification_table = slopes_lin_amplification;
save('slopes_threshold_sigma_0_large_full','threshold_lin_amp_table','threshold_zero_force_table','slopes_threshold_lin_amp_table','slopes_threshold_zero_force_table','slopes_lin_amplification_table','q_table');








disp('done')



