function thresh_out  = soft_threshold(val, threshold)

thresh_out = sign(val) .* max( (abs(val)- threshold) ,0 );