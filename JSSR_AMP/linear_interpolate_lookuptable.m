function y_value = linear_interpolate_lookuptable(x_value,x_table,y_table)

%function y_value = linear_interpolate_lookuptable(x_value,x_table,y_table)

ind_smaller = find(x_table <= x_value);
if( isempty(ind_smaller))
    y_value = y_table(1);
    return
end
size_ind_smaller = size(ind_smaller);

ind_greater = find(x_table >= x_value);
if(isempty(ind_greater))
    table_size = size(y_table);
    y_value = y_table(table_size(2));
    return
end
size_ind_greater = size(ind_greater);

x_value_lower_ind = x_table(ind_smaller(size_ind_smaller(2)));
x_value_higher_ind = x_table(ind_greater(1));

y_value_lower_ind = y_table(ind_smaller(size_ind_smaller(2)));
y_value_higher_ind = y_table(ind_greater(1));

if(x_value_lower_ind == x_value_higher_ind)
    y_value = y_table(find(x_table == x_value));
else
    y_value = y_value_lower_ind + ( y_value_higher_ind - y_value_lower_ind)/(x_value_higher_ind - x_value_lower_ind)*(x_value-x_value_lower_ind);
end 
    

