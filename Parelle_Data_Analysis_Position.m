num_of_data = OutputNumber * core_number;
px = 0;
py = 0;
final_mean = zeros([1,TempSlice]);
final_svar = zeros([1,TempSlice]);

for time_index = 1:1:TempSlice
    data_batch = [];
    for zjy_index = 1:1:core_number
        for output_index = 1:1:OutputNumber
            for site_index = 1:1:NumOfVertexs
                data_batch(end+1) = propa_green_group(site_index,site_index,time_index,output_index,zjy_index);
            end
        end
    end
    %final_mean(time_index) = mean(data_batch);
    %final_svar(time_index) = sqrt(var(data_batch));
    [final_mean(time_index),final_svar(time_index)] = Bin(data_batch);
end
%errorbar(final_mean,final_svar) %%Before Log
Tau_range = D_Tau:D_Tau:Beta;
log_mean = log(final_mean);
log_svar = final_svar./final_mean;
errorbar(Tau_range,final_mean,final_svar,'r');
title(['L=',num2str(NumInEdge),'   ','p_x=',num2str(px),'   ','p_y=',num2str(py),'   \beta=',num2str(Beta),'   ','Uene=',num2str(Uene),'   ','\Delta_\tau=',num2str(D_Tau)])
xlabel("\tau");
ylabel("G(\tau)")

%% Then we calculate the theoretical result for position space
Miu = Uene/2;
plot_sum = 0.0;
for px_index = 1:1:NumInEdge
    for py_index = 1:1:NumInEdge
        delta_p = 2*pi/NumInEdge;
        px = px_index*delta_p-pi;
        py = py_index*delta_p-pi;
        e_k = - 2*T_hop*(cos(px)+cos(py));
        pre_factor = 1.0-1.0/(exp(Beta*(e_k-Miu))+1);
        plot_theo = pre_factor * exp(-e_k*Tau_range);
        plot_sum = plot_sum + plot_theo/NumOfVertexs;
    end
end
%plot(Tau_range,plot_sum);
%hold on
%plot(Tau_range,final_mean,'r*')


%% And here the below is a virtual analytic continuation
delta = 0.4;
plot_sum = 0.0;
ene_range = -20:0.1:20;
for px_index = 1:1:NumInEdge
    for py_index = 1:1:NumInEdge
        delta_p = 2*pi/NumInEdge;
        px = px_index*delta_p-pi;
        py = py_index*delta_p-pi;
        e_k = - 2*T_hop*(cos(px)+cos(py));
        spec_single = 1.0/pi * (delta./((ene_range-e_k).^2 + delta.^2));
        plot_sum = plot_sum + spec_single/(NumInEdge)^2;
    end
end
plot(ene_range,plot_sum,'r*')
hold on
plot(ene_range,plot_sum);
    