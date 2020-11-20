N = fix(length(mea_result_auxi)*0.1);
sample = 2*mea_result_auxi(1:N);
%N = 1e4;
%sample = rand([1,N])+1;
%%%%%%% JackKnife %%%%%%%%
pseudo_value = zeros([1,N]);
ori_mean = mean(sample);
for index = 1:1:N
    if mod(index,1000)==0
        disp(index/N);
    end
    temp_jack = sample;
    temp_sum = sum(sample)-sample(index);
    pseudo_mean = temp_sum/(N-1);
    %pseudo_value_ans = N*ori_mean - (N-1)*pseudo_mean;
    pseudo_value_ans = pseudo_mean;
    pseudo_value(index) = pseudo_value_ans;
end
jack_mean = mean(pseudo_value);
jack_var = sum((pseudo_value - jack_mean).^2);
jack_svar = sqrt(jack_var);
%%%%%%%%%%%%%% Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%
boot_mean_list = zeros([1,N]);
for index = 1:1:N
    if mod(index,100)==0
        disp(index/N);
    end
    boot_temp = randi([1,N],[1,N]);
    for temp_index = 1:1:N
        boot_temp(temp_index) = sample(boot_temp(temp_index));
    end
    boot_mean_list(index) = mean(boot_temp);
end
boot_svar = sqrt(mean(boot_mean_list.^2) - mean(boot_mean_list)^2);
        



%%%%%%%%%%%% Then we play we bin method %%%%%%%%%%%%%%%%%
bin_step = 1;
bin_svar_data = [];
bin_size_data = [];
for bin_size = bin_step:bin_step:N/10
    if mod(N,bin_size)==0
        bin_size_data(end+1) = bin_size;
        bin_num = N/bin_size;
        bin_data = zeros([1,bin_num]);
        for index = 1:1:bin_num
           bin_data(index) = mean(sample(bin_size*(index-1)+1:bin_size*index)); 
        end
        bin_mean = mean(bin_data);
        bin_var = sum((bin_data - bin_mean).^2)/(bin_num);
        bin_svar = sqrt( bin_var );
        bin_svar_data(end+1) = bin_svar;
    end
end
plot(bin_size_data,bin_svar_data);
hold on;
plot(bin_size_data,bin_svar_data,'r*');
hold on
xlabel("BinSize");
ylabel("Standard Variance");
title("StdVar for Double Occupancy")
plot(bin_size_data,0.0*bin_size_data + jack_svar,'m');
hold on
plot(bin_size_data,0.0*bin_size_data + boot_svar,'r');
legend(["Bin Method","Bin Method","Jackknife Method","BootStrap"]);