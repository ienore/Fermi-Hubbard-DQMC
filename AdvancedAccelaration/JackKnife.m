function [jack_mean,jack_stdvar] = JackKnife(sample,N)
pseudo_value = zeros([1,N]);
for index = 1:1:N
    temp_sum = sum(sample)-sample(index);
    pseudo_mean = temp_sum/(N-1);
    %pseudo_value_ans = N*ori_mean - (N-1)*pseudo_mean;
    pseudo_value_ans = pseudo_mean;
    pseudo_value(index) = pseudo_value_ans;
end
jack_mean = mean(pseudo_value);
jack_var = sum((pseudo_value - jack_mean).^2);
jack_stdvar = sqrt(jack_var);
end

