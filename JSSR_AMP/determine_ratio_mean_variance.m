%function [ratio_maxactual_maxapprox q_table] = determine_ratio_maxactual_maxapprox(sigma_0,sigma_1)

clc
clear all;

large_number =10000;

sigma_1=10;
sigma_0=1;

mean_1=0;
mean_0=0;

variance_1 = sigma_1*sigma_1;
variance_0 = sigma_0*sigma_0;

step_size_x =0.01;
x= -50:step_size_x:50;


step_size_q =0.05;
init_val_q = 0.1;
range_q = 0.95;

num_of_steps_q = (range_q - init_val_q)/step_size_q;

for i=1:num_of_steps_q % for each value of q
    
q1= init_val_q + step_size_q*(i-1);
disp(q1);

q(i)= q1;

q0=1-q1;
    
range_upper_mean_prod =30;
mean_prod_step_size = 0.25;
num_of_steps_mean_prod = range_upper_mean_prod/mean_prod_step_size;
%num_of_steps_mean_prod = 5;

for k=1:1:num_of_steps_mean_prod % for each value of mean
    
    disp('mean = '); 
    mean_prod =mean_prod_step_size*k;
    mean_prod_x_axis(k) = mean_prod;
    disp(mean_prod);
    
range_upper_variance = 45;
variance_step_size = 0.01;
num_of_steps_variance = range_upper_variance/variance_step_size;
%num_of_steps_variance = 40;

for j=1:1:num_of_steps_variance % for each variance
    
variance_prod =  j*variance_step_size;
var(j) = variance_prod;
sigma_prod = sqrt(variance_prod);






[weight_0_new  mean_final_0 variance_final_0] = weight_determination(variance_0,mean_0,variance_prod,mean_prod,q0);
[weight_1_new  mean_final_1 variance_final_1] = weight_determination(variance_1,mean_1,variance_prod,mean_prod,q1);
           
            small_index_weight_0 =  find(weight_0_new<=1e-30);
            weight_0_new(small_index_weight_0) = 0;
            weight_1_new(small_index_weight_0)= 1;
           
            weight_0_final = weight_0_new ./ (weight_0_new + weight_1_new);
            weight_1_final = weight_1_new ./ (weight_0_new + weight_1_new);
            
            weight_0(j,k,i)=weight_0_final;
            weight_1(j,k,i)=weight_1_final;
 
            
            approx_variance_including_weights(j) = weight_0_final.*variance_final_0 + weight_1_final.*variance_final_1 +...
                weight_0_final.*weight_1_final.*(mean_final_0 - mean_final_1).^2;
            
            

                
         pdf_max_value_1 = weight_0_new * normpdf(mean_final_0,mean_final_0,sqrt(variance_final_0)) + ...
                           weight_1_new * normpdf(mean_final_0,mean_final_1,sqrt(variance_final_1)); 
         pdf_max_value_2 = weight_1_new * normpdf(mean_final_1,mean_final_1,sqrt(variance_final_1)) + ...
                           weight_0_new * normpdf(mean_final_1,mean_final_0,sqrt(variance_final_0)); 
                     
         if (pdf_max_value_1 > pdf_max_value_2)
             mean_actual(j)= mean_final_0;
         else
             mean_actual(j)= mean_final_1;
         end  

         approx_variance(j) = (q0*variance_0 + q1*variance_1)* variance_prod /(q0*variance_0 + q1*variance_1 + variance_prod );
         mean_approx(j) = approx_variance(j)* mean_prod / variance_prod;
         
         
         
         
         
         
         
                   
end

ratio_test_mean(:,i,k) = mean_actual./mean_approx;
ratio_mean(:,k,i) = mean_actual./mean_approx;
ratio_var(:,k,i) = approx_variance_including_weights./approx_variance;
ratio_test_var(:,i,k) = approx_variance_including_weights./approx_variance;


end

end
q_table =q;

save('ratio_maxactual_max_approx_large_full','ratio_mean','ratio_var','q_table','range_upper_variance','variance_step_size','range_upper_mean_prod','mean_prod_step_size');
