px_index = 2;
py_index = 1;
N_epoch = 1e4;
eps = 1e-3;
step = 2e-0;
E_min = -10;
E_max = 10;
dE = 0.2;
alpha = 0.05;
E_range = E_min:dE:E_max;

G_k = zeros([TempSlice,1]);
G_k = propa_mean';
Gauge = zeros([TempSlice,length(E_range)]);



for time_index = 1:1:TempSlice
    for E_index = 1:1:length(E_range)
        tau = (time_index-1)*D_Tau;
        omega = E_range(E_index);
        Gauge(time_index,E_index) = dE*exp(-omega*(tau-Beta/2))/(2*cosh(Beta*omega/2));
    end
end

A = zeros([length(E_range),1]);



for epoch_index = 1:1:N_epoch
    if mod(epoch_index,10)==0
        fprintf(" Epoch_ratio = %f, diff = %f\n",epoch_index/N_epoch,diff)
    end
    for A_index = 1:1:length(E_range)
        A_new = A;
        A_new(A_index) = A_new(A_index) + step*(rand-0.5);
        diff = norm(Gauge*A - G_k);
        value_old = sum( (Gauge*A - G_k).^2 ) + alpha * dE* sum( (A(1:(length(A)-1))-A(2:length(A))).^2 );
        value_new = sum((Gauge*A_new - G_k).^2) + alpha * dE* sum( (A_new(1:(length(A)-1))-A(2:length(A))).^2 );
        if value_new<value_old && A_new(A_index)>=0;
            A = A_new;
        end
    end
end

% delta_p = 2*pi/(NumInEdge);
% px = (px_index-NumInEdge) * delta_p + pi ;
% py = (py_index-NumInEdge) * delta_p + pi;
px = 0;
py = pi;
 
plot(E_range,A,'r*')
hold on
plot(E_range,A)
xlabel("\omega");
ylabel("A(\omega)")
title(['Specturm  Uene = ',num2str(Uene),'  ','Beta = ',num2str(Beta),'   ','D_Tau = ',num2str(D_Tau), ...
    'px = ',num2str(px),'  ','py = ',num2str(py)])

