N_epoch = 1e5;
eps = 1e-3;
step = 1e-2;
E_min = -10;
E_max = 10;
dE = 0.05;
E_range = E_min:dE:E_max;

CHANGE_BOUND = 2e3;
alpha = 0.01;
alpha_decay_ratio = 0.2;
%G_k = zeros([TempSlice,1]);
Gauge = zeros([TempSlice,length(E_range)]);
count = 0;
diff_in_Xi = 1.0;
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
for time_index = 1:1:TempSlice+1
    for E_index = 1:1:length(E_range)
        tau = (time_index-1)*D_Tau;
        omega = E_range(E_index);
        Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
    end
end
plot_mean_ori = final_mean;
lenG = length(plot_mean_ori);
G_k = zeros([lenG+1,1]);
G_k(2:lenG+1) = plot_mean_ori(1:lenG);
G_k(1) = plot_mean_ori(lenG);

norm_Gauge = norm(Gauge);
A = rand([length(E_range),1]);
change_sign = 1;
diff_in_Xi = 1.0;
diff_in_Ent = 0.0;
count_not_change = 0;
epoch_index = 0;
while (epoch_index < N_epoch || count_not_change < CHANGE_BOUND-10)
    epoch_index = epoch_index + 1;
    if mod(epoch_index,N_epoch/500)==0 || epoch_index == 1
        temp_vec = Gauge*2*A - 2*G_k;
        Err = max(abs(Gauge*A-G_k));
        step = Err/norm_Gauge;
        fprintf("Epoch Ratio = %f,  diff = %e, alpha = %e\n",epoch_index/N_epoch,max(abs((Gauge*A-G_k)./G_k)),alpha);
        if count_not_change > CHANGE_BOUND
            alpha = alpha_decay_ratio * alpha;
        end
    end
    for A_index = 1:1:length(E_range)
        dA = step*(rand-0.5);
        A_pre = A(A_index);
        diff_in_Xi = (Gauge(:,A_index)*dA)'*(temp_vec +  Gauge(:,A_index)*dA);
        Ent_old = dE*A_pre*log(A_pre/(1.0/(E_max-E_min)));
        Ent_new = dE*(A_pre+dA)*log((A_pre+dA)/(1.0/(E_max-E_min)));
        diff_in_Ent = alpha*(Ent_new-Ent_old);
        if diff_in_Xi+diff_in_Ent < 0 && (A_pre + dA)>0
            A(A_index) = A_pre+dA;
            temp_vec = temp_vec + 2*Gauge(:,A_index)*dA;
            count_not_change = 0;
        else
            count_not_change = count_not_change +1;
        end
    end
end
plot(E_range,A,'r*')
hold on
plot(E_range,A,'r')
diff_G = (Gauge*A - G_k)./G_k;
xlabel("\omega")
ylabel("A(\omega)")
title(['L=',num2str(NumInEdge),'   ','  \beta=',num2str(Beta),'   ','Uene=',num2str(Uene),'   ','\Delta_\tau=',num2str(D_Tau)])

%plot(diff_G,'b')