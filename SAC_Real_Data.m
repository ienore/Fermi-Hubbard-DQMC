
N_epoch = 1e4;
eps = 1e-3;
step = 1e-2;
E_min = -5;
E_max = 20;
dE = 0.5;
E_range = E_min:dE:E_max;
alpha = 0.1;
G_k = zeros([length(plot_mean_ori)+1,1]);
G_k(2:length(plot_mean_ori)+1) = plot_mean_ori';
G_k(1) = G_k(length(plot_mean_ori)+1);
Gauge = zeros([TempSlice+1,length(E_range)]);
for time_index = 1:1:TempSlice+1
    for E_index = 1:1:length(E_range)
        tau = (time_index)*D_Tau;
        omega = E_range(E_index);
        Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
    end
end


A = zeros([length(E_range),1]);
for epoch_index = 1:1:N_epoch
    if mod(epoch_index,N_epoch/1000)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        A_new = A_new/(sum(A_new)*dE);
        diff = norm(Gauge*A - G_k);
        value_old = sum( (Gauge*A - G_k).^2 ) + alpha * dE* sum( (A(1:(length(A)-1))-A(2:length(A))).^2 );
        value_new = sum((Gauge*A_new - G_k).^2) + alpha * dE* sum( (A_new(1:(length(A)-1))-A(2:length(A))).^2 );
        if value_new<value_old && A_new(A_index) >= 0
            A = A_new;
            %A = A/(sum(A)*dE);
        end
    end
end
 
%plot(E_range,A,'r*')
%hold on
%plot(E_range,A)
diff_G = (Gauge*A - G_k)./G_k;
%plot(diff_G,'b')
%plot(E_range,A_ori,'r');
%hold on
plot(E_range,A,'b')



