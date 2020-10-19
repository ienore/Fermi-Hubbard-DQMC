num_of_data = OutputNumber * core_number;
px = 0;
py = pi;
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
plot_mean(2:TempSlice+1) = final_mean;
plot_svar(2:TempSlice+1) = final_svar;
plot_mean(1) = final_mean(TempSlice);
plot_svar(1) = final_svar(TempSlice);
errorbar(plot_mean,plot_svar,'r');
title(['L=',num2str(NumInEdge),'   ','p_x=',num2str(px),'   ','p_y=',num2str(py),'   \beta=',num2str(Beta),'   ','Uene=',num2str(Uene),'   ','\Delta_\tau=',num2str(D_Tau)])
    
    