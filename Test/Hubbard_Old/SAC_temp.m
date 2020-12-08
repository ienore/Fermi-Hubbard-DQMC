px_index = 1;
py_index = 1;
N_epoch = 1e3;
eps = 1e-3;
step = 2e-1;
E_min = -10;
E_max = 10;
dE = 1;
E_range = E_min:dE:E_max;

G_k = zeros([TempSlice,1]);
G_k(2:TempSlice) = mea_ML_data_p(px_index,py_index,1:TempSlice-1);
G_k(1) = mea_ML_data_p(px_index,py_index,1);
Gauge = zeros([TempSlice,length(E_range)]);


A = zeros([length(E_range),1]);
A(7) = 1;
A(13) = 1;
A_init = A;


for time_index = 1:1:TempSlice
    for E_index = 1:1:length(E_range)
        tau = (time_index-1)*D_Tau;
        omega = E_range(E_index);
        Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
    end
end
G_k = Gauge*A;


A = rand([length(E_range),1]);
for epoch_index = 1:1:N_epoch
    if mod(epoch_index,10)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        diff = norm(Gauge*A_new - G_k);
        if diff < norm(Gauge*A - G_k) && A_new(A_index)>0;
            A = A_new;
        end
    end
end
G_diff = Gauge*A - G_k;
fprintf("Diff = %f\n",diff);
plot(G_diff);
% plot(E_range,A,'r*')
% hold on
% plot(E_range,A)
% xlabel("\omega");
% ylabel("A(\omega)")
% title(['Specturm  Uene = ',num2str(Uene),'  ','Beta = ',num2str(Beta),'   ','D_Tau = ',num2str(D_Tau)])

G_init = Gauge * A_init;
G_k = Gauge * A;
plot(G_init - G_k);
plot(A_init-A);
