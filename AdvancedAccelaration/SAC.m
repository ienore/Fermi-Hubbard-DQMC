
N_epoch = 1e4;
eps = 1e-3;
step = 2e-2;
E_min = -5;
E_max = 5;
dE = 0.25;
E_range = E_min:dE:E_max;
A_mat = zeros([length(E_range),NumInEdge,NumInEdge]);
count = 0;
for px_index = 1:1:NumInEdge
    for py_index = 1:1:NumInEdge
        count = count + 1;
        G_k = zeros([TempSlice,1]);
        G_k(2:TempSlice) = mea_ML_data_p(px_index,py_index,1:TempSlice-1);
        G_k(1) = mea_ML_data_p(px_index,py_index,1);
        Gauge = zeros([TempSlice,length(E_range)]);

        for time_index = 1:1:TempSlice
            for E_index = 1:1:length(E_range)
                tau = (time_index-1)*D_Tau;
                omega = E_range(E_index);
                Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
            end
        end
        A = rand([length(E_range),1]);
        for epoch_index = 1:1:N_epoch
            if mod(epoch_index,1000)==0
                fprintf("Count = %d, Epoch_ratio = %f, diff = %f\n",count,epoch_index/N_epoch,diff)
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
        A_mat(:,px_index,py_index) = A(:);
    end
end
A_final = zeros([length(E_range),1]);
for px_index = 1:1:NumInEdge
    for py_index = 1:1:NumInEdge
        A_final = A_final + A_mat(:,px_index,py_index)/NumOfVertexs;
    end
end
plot(E_range,A_final,'r*')
hold on
plot(E_range,A_final)
xlabel("\omega");
ylabel("A(\omega)")
title(['Specturm  Uene = ',num2str(Uene),'  ','Beta = ',num2str(Beta),'   ','D_Tau = ',num2str(D_Tau)])