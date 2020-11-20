N = 1e5;
T = 0.1;
x_old = 0;
RA = 50;%rou_max = N/RA
sample = zeros([1,N]);
step = 200;
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
mean_sample = mean(sample);
rou_list = zeros([1,fix(N/RA)]);
bin_range = 1:1:fix(N/RA);
for index_rou = 1:1:fix(N/RA)
    if mod(index_rou,fix(N/(RA*100)))==0
        disp(index_rou/(N/RA));
    end
    rou_sum = mean((sample(1:(N-index_rou)).* sample((index_rou+1):N)));
    rou_list(index_rou) = rou_sum;
end
rou_list = rou_list/rou_list(1);
rou_list_plot = [];
for index = 1:1:length(rou_list)
    if rou_list(index) < 0.2
        break;
    end
end
rou_list_plot = rou_list(1:index);   
rou_range = 1:1:length(rou_list_plot);

%plot(rou_range,log(rou_list_plot),'r*');
%hold on
%plot(rou_range,log(rou_list_plot),'b');
tau_int = sum(rou_list);
%%%%%%%%%%%%%%%%%%%% The Bin Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin_step = 1e1;
N_bin = 1e4;
bin_err_data = zeros([1,N_bin]);
bin_err_xrange = (1:1:N_bin)*bin_step;
jack_err_data = zeros([1,N_bin]);
bin_result = zeros([1,N_bin]);
for index = 1:1:N_bin
    if mod(index,N_bin/100)==0
        disp(index/N_bin);
    end
    %disp(index);
    bin_size = bin_err_xrange(index);
    max_bin_num = fix(N/bin_size);
    bin_temp = zeros([1,max_bin_num]);
    for temp_index = 1:1:max_bin_num
       temp = sample((temp_index-1)*bin_size+1:temp_index*bin_size);
       bin_temp(temp_index) =  mean(temp);
    end
    bin_result(index) = mean((temp - mean(temp)).^2);
    bin_err_data(index) = sqrt((sum(bin_temp.^2)/max_bin_num - (sum(bin_temp)/max_bin_num)^2)/(max_bin_num-1));
end

for index = 1:1:length(bin_err_data)
   tau = bin_err_data(index)^2/(2*bin_err_data(1)^2);
   plot_new(index) = bin_err_xrange(index)/tau;
end


plot(bin_err_xrange,bin_result/T,'r')
