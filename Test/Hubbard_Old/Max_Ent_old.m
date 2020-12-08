
TempSlice = 32;
Beta = 2.0;
D_Tau = Beta / TempSlice;
N_epoch = 1e4;
eps = 1e-3;
step = 1e-2;
E_min = -15;
E_max = 15;
dE = 0.2;
E_range = E_min:dE:E_max;
alpha = 0.01;
G_k = zeros([TempSlice,1]);
Gauge = zeros([TempSlice,length(E_range)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros([length(E_range),1]);
for index = 1:1:length(E_range)
    omega = E_range(index);
    A(index) = exp(-0.2*(omega+5)^2)+exp(-0.2*(omega-5)^2);
    %A(index) = exp(-0.2*omega^2);
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
%G_k = propa_mean';
A = rand([length(E_range),1]);
for epoch_index = 1:1:N_epoch
    if mod(epoch_index,10)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        %A_new = A_new/(sum(A_new)*dE);
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