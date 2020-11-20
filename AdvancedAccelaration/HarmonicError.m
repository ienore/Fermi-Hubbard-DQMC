N = 1e5;
T = 0.1;
x_old = 0;
step = 200;%增大这个，整体的关联长度就会上升
sample = zeros([1,N]);
count = 1;
for iter = 1:1:N
    x_new = x_old + (rand-0.5)*step;
    if 0.5*x_new^2 < 0.5*x_old^2
        x_old = x_new;
    else
        if rand < exp(- 1/T *(0.5*x_new^2 - 0.5*x_old^2))
            x_old = x_new;
        end
    end
    sample(count) = x_old;
    count = count + 1;
end
sample = sample;



%%%%%%% JackKnife %%%%%%%%
% sample_jack = sample;
% N_jack = length(sample_jack);
% jack_data = zeros([1,N_jack]);
% all_sum = sum(sample_jack);
% for jack_index = 1:1:N_jack
%     jack_data(jack_index)= (all_sum - sample_jack(jack_index))/(N_jack-1);
% end
% jack_err_mean = mean(jack_data);
% jack_diff = jack_data - jack_err_mean;
% jack_err = sqrt((N_jack-1)/N_jack*sum(jack_diff.^2));
% 
% 

% %%%%%%%%%%%%%% Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boot_mean_list = zeros([1,N]);
% for index = 1:1:N
%     if mod(index,100)==0
%         disp(index/N);
%     end
%     boot_temp = randi([1,N],[1,N]);
%     for temp_index = 1:1:N
%         boot_temp(temp_index) = sample(boot_temp(temp_index));
%     end
%     boot_mean_list(index) = mean(boot_temp);
% end
% boot_svar = sqrt(mean(boot_mean_list.^2) - mean(boot_mean_list)^2);
%         



%%%%%%%%%%%% Correlation Function %%%%%%%%%%%%%%%%%

% plot_step = 1e0;
% N_corr = 2000;
% count = 1;
% correlation_data = zeros([1,N_corr]);
% dis_range = (1:1:N_corr) * plot_step;
% for plot_index = 1:1:N_corr
%     dis = plot_index*plot_step;
%     temp_1 = sample(1:N-dis);
%     temp_2 = sample(dis+1:N);
%     correlation_data(count) = (temp_1*transpose(temp_2))/length(temp_1);
%     count = count + 1;
% end
% plot(correlation_data,'r')
% xlabel("Distance Between");
% ylabel("Correlation Function for x");
% title(['Correlation Function','   ','Step=',num2str(step)]);



%%%%%%%%%%%%%%%%%%%% The Bin Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin_step = 1e0;
N_bin = 2e3;
bin_err_data = zeros([1,N_bin]);
bin_err_xrange = (1:1:N_bin)*bin_step;
jack_err_data = zeros([1,N_bin]);
for index = 1:1:N_bin
    if mod(index,N_bin/100)==0
        disp(index/N_bin);
    end
    %disp(index);
    bin_size = index * bin_step;
    max_bin_num = fix(N/bin_size);
    bin_temp = zeros([1,max_bin_num]);
    for temp_index = 1:1:max_bin_num
       temp = sample((temp_index-1)*bin_size+1:temp_index*bin_size);
       bin_temp(temp_index) =  mean(temp);
    end
    bin_err_data(index) = sqrt((sum(bin_temp.^2)/max_bin_num - (sum(bin_temp)/max_bin_num)^2)/(max_bin_num-1));
    sample_jack = bin_temp;
    N_jack = length(sample_jack);
    jack_data = zeros([1,N_jack]);
    all_sum = sum(sample_jack);
    for jack_index = 1:1:N_jack
        jack_data(jack_index)= (all_sum - sample_jack(jack_index))/(N_jack-1);
    end
    jack_err_mean = mean(jack_data);
    jack_diff = jack_data - jack_err_mean;
    jack_err_data(index) = sqrt((N_jack-1)/N_jack*sum(jack_diff.^2));
end
plot(bin_err_xrange,bin_err_data,'r')
hold on
%plot(bin_err_xrange,bin_err_xrange*0+jack_svar);
xlabel("Distance Between");
ylabel("Err for x");
title(['Err','   ','Step=',num2str(step)]);
%plot(bin_err_xrange,err_diff,'r')
%legend("BinMethod","JackKnife")
legend("Err for Bin Methods");

% 
% plot(bin_size_data,bin_svar_data);
% hold on;
% plot(bin_size_data,bin_svar_data,'r*');
% hold on
% xlabel("BinSize");
% ylabel("Standard Variance");
% title("StdVar for Double Occupancy")
% plot(bin_size_data,0.0*bin_size_data + jack_svar,'m');
% hold on
% plot(bin_size_data,0.0*bin_size_data + boot_svar,'r');
% legend(["Bin Method","Bin Method","Jackknife Method","BootStrap"]);