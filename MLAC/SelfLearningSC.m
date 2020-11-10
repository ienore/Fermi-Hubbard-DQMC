
TempSlice = 16;
Beta = 2.0;
D_Tau = Beta / TempSlice;
N_epoch = 1e4;
N_warm = 1e4;
eps = 1e-3;
step = 1e-2;
E_min = -20;
E_max = 20;
dE = 0.2;
E_range = E_min:dE:E_max;
alpha = 0.02;
G_k = zeros([TempSlice,1]);
Gauge = zeros([TempSlice,length(E_range)]);
Theta = 1e5;%Annealing Factor
%%%%%%%%%%%%%%%%%%%%% Set Up A_origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros([length(E_range),1]);
for index = 1:1:length(E_range)
    omega = E_range(index);
    %A(index) = exp(-0.2*(omega+5)^2)+exp(-0.2*(omega-5)^2);
    A(index) = exp(-0.2*omega^2);
end
A = A/(sum(A)*dE);
A_ori = A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time_index = 1:1:TempSlice
    for E_index = 1:1:length(E_range)
        tau = (time_index-1)*D_Tau;
        omega = E_range(E_index);
        Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
    end
end
G_k = Gauge*A;
A = rand([length(E_range),1]);
A = A/(sum(A)*dE);

%%%%%%%%%%%%%%%%%%%%%%%% Calculation Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_collect = zeros([length(E_range),N_epoch]);
A = A/(sum(A)*dE);
for epoch_index = 1:1:N_warm
    if mod(epoch_index,100)==0
        fprintf(" Warm_ratio = %f, diff = %f\n",epoch_index/N_warm,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        A_new = A_new/(sum(A_new)*dE);
        diff = norm(Gauge*A - G_k);
        value_old = sum( (Gauge*A - G_k).^2 );
        value_new = sum((Gauge*A_new - G_k).^2);
        step = sqrt(diff);
        if rand < exp(Theta*(value_old-value_new)) && A_new(A_index) >= 0;
            A = A_new;
        end
    end
end
for epoch_index = 1:1:N_epoch
    if mod(epoch_index,100)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        A_new = A_new/(sum(A_new)*dE);
        diff = norm(Gauge*A - G_k);
        value_old = sum( (Gauge*A - G_k).^2 );
        value_new = sum((Gauge*A_new - G_k).^2);
        step = sqrt(diff);
        if rand < exp(Theta*(value_old-value_new)) && A_new(A_index) >= 0;
            A = A_new;
        end
    end
    A_collect(:,epoch_index) = A;
end
A_final = zeros([length(E_range),1]);
for index = 1:1:length(E_range);
   A_final(index) = mean(A_collect(index,:));  
end
diff_G = (Gauge*A_final - G_k)./G_k;
%plot(diff_G,'b')
plot(E_range,A_ori,'r');
hold on
plot(E_range,A_final,'b')

