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
errorbar(final_mean,final_svar)
    
    