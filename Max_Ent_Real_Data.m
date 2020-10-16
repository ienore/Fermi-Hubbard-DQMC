
N_epoch = 1e4;
eps = 1e-3;
step = 1e-2;
E_min = -10;
E_max = 10;
dE = 0.1;
E_range = E_min:dE:E_max;
alpha = 0.01;

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


A = rand([length(E_range),1]);
for epoch_index = 1:1:N_epoch
    if mod(epoch_index,10)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        A_new = A_new/(sum(A_new)*dE);
        diff = norm(Gauge*A - G_k);
        step = sqrt(diff);
        G_new = Gauge*A_new;
        G = Gauge*A;
        value_old = sum( (G - G_k).^2 ) + alpha * dE* sum(A.*log(A/(1.0/(E_max-E_min))));
        value_new = sum((G_new - G_k).^2) +alpha * dE* sum(A_new.*log(A_new/(1.0/(E_max-E_min))));
        if value_new<value_old && A_new(A_index) >= 0;
            A = A_new;
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

