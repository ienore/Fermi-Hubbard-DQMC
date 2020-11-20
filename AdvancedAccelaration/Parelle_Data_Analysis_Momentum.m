num_of_data = OutputNumber * core_number;
px = pi/2;
py = pi/2;
final_mean = zeros([1,TempSlice]);
final_svar = zeros([1,TempSlice]);

for time_index = 1:1:TempSlice
    data_batch = [];
    for zjy_index = 1:1:core_number
        for output_index = 1:1:OutputNumber
            temp_green_sum = 0.0;
            for x_index = 1:1:NumOfVertexs
                for y_index = 1:1:NumOfVertexs
                    [x_a,y_a] = IndexToCoor(x_index,NumInEdge);
                    [x_b,y_b] = IndexToCoor(y_index,NumInEdge);
                    temp_green_sum = temp_green_sum + exp(1i * (px*(x_b-x_a)+py*(y_b-y_a))) ...
                    *propa_green_group(x_index,y_index,time_index,output_index,zjy_index);
                    %disp(propa_green_group(x_index,y_index,time_index,output_index,zjy_index));
                end
            end
            data_batch(end+1) = real(temp_green_sum)/NumOfVertexs;
        end
    end
    %final_mean(time_index) = mean(data_batch);
    %final_svar(time_index) = sqrt(var(data_batch));
    [final_mean(time_index),final_svar(time_index)] = Bin(data_batch);
end
plot_mean = zeros([1,TempSlice+1]);
plot_svar = zeros([1,TempSlice+1]);
%% Simply Get G without any modification
plot_mean = final_mean;
plot_svar = final_svar;
%% Get G(\beta) to be G(0)
% plot_mean(2:TempSlice+1) = final_mean;
% plot_svar(2:TempSlice+1) = final_svar;
% plot_mean(1) = final_mean(TempSlice);
% plot_svar(1) = final_svar(TempSlice);
Tau_range = D_Tau:D_Tau:Beta;
%errorbar(Tau_range,plot_mean,plot_svar,'r');
errorbar(Tau_range,plot_mean,plot_svar);
%title(['L=',num2str(NumInEdge),'   ','p_x=',num2str(px),'   ','p_y=',num2str(py),'   \beta=',num2str(Beta),'   ','Uene=',num2str(Uene),'   ','\Delta_\tau=',num2str(D_Tau)])
xlabel("\tau");
ylabel("G(k,\tau)")
hold on
%% Then we calculate the theoretical result
Miu = Uene/2;
e_k = - 2*T_hop*(cos(px)+cos(py));
pre_factor = 1.0-1.0/(exp(Beta*(e_k-Miu))+1);
plot_theo = pre_factor * exp(-e_k*Tau_range);
plot(Tau_range,plot_theo,'blue');
hold on
errorbar(Tau_range,(plot_mean),plot_svar,'r*');


%% Formal Plot
errorbar(Tau_range,(plot_mean),plot_svar,'r');
title(['L=',num2str(NumInEdge),'   ','p_x=',num2str(px),'   ','p_y=',num2str(py),'   \beta=',num2str(Beta),'   ','Uene=',num2str(Uene),'   ','\Delta_\tau=',num2str(D_Tau)])
xlabel("\tau");
ylabel("G(k,\tau)")

    