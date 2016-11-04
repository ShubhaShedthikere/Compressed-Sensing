function [weight_new mean_final variance_final] = weight_determination(rep_var,rep_mean,ret_prod_incoming_messages_variance,ret_prod_incoming_message_mean,rep_weight)

variance_final = rep_var.*ret_prod_incoming_messages_variance ./ (rep_var + ret_prod_incoming_messages_variance);


mean_final =  variance_final.*((ret_prod_incoming_message_mean./ret_prod_incoming_messages_variance));

var_sum_sqrt_inv = 1./sqrt(rep_var + ret_prod_incoming_messages_variance);

exponent = exp(0.5 *min((mean_final.*mean_final ./ variance_final),500)).* ...
           exp( - 0.5*( min(rep_mean.*rep_mean ./rep_var,0) + ret_prod_incoming_message_mean.*ret_prod_incoming_message_mean./ret_prod_incoming_messages_variance));

weight_new = (rep_weight/sqrt(2*pi)).* var_sum_sqrt_inv .* exponent;
%weight_new = rep_weight;


