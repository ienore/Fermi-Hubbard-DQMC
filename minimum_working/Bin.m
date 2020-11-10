function [bin_mean,bin_stdvar] = Bin(sample)
    N = length(sample);
    bin_num = 5;
    max_bin_num = bin_num;
    bin_size = fix(N/bin_num);
    bin_data = zeros([1,bin_num]);
    bin_mean = mean(sample);
    for index = 1:1:bin_num
        temp = sample((index-1)*bin_size+1:index*bin_size);
        bin_data(index) = mean(temp);
    end
    bin_temp = bin_data - bin_mean;
    bin_stdvar = sqrt((sum(bin_temp.^2)/max_bin_num - (sum(bin_temp)/max_bin_num)^2)/(max_bin_num-1)); 
    %disp(bin_stdvar);
end

