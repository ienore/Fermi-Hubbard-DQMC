%sample = rand([1,1e4]);
sample = sample - mean(sample);

T = 0.1;
x_old = 0;
step = 200;%增大这个，整体的关联长度就会上升
count = 1;
N = length(sample);


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
N_bin = 1e2;
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
end
plot(bin_err_xrange,bin_err_data,'r')
hold on
%plot(bin_err_xrange,bin_err_xrange*0+jack_svar);
xlabel("Distance Between");
ylabel("Err for x");
title(['Err','   ','Step=',num2str(step)]);