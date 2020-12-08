green_propa_final = zeros([NumOfVertexs,NumOfVertexs,TempSlice]);
green_propa_final_svar = zeros([NumOfVertexs,NumOfVertexs,TempSlice]);
for x_index = 1:1:NumOfVertexs
    for y_index = 1:1:NumOfVertexs
        for t_index = 1:1:TempSlice
            buff_pre_calc = reshape(propa_green_group(x_index,y_index,t_index,:,:),[1,OutputNumber*NumOfCores]);
            green_propa_final(x_index,y_index,t_index) = mean(buff_pre_calc);
            green_propa_final_svar(x_index,y_index,t_index) = sqrt(var(buff_pre_calc));
        end
    end
end

sum_diag = zeros([1,TempSlice]);
sum_diag_svar = zeros([1,TempSlice]);
for t_index = 1:1:TempSlice
    temp = zeros([1,NumOfVertexs]);
    temp_svar = zeros([1,NumOfVertexs]);
    for site_index = 1:1:NumOfVertexs
        temp(site_index) = green_propa_final(site_index,site_index,t_index);
        temp_svar(site_index) = green_propa_final_svar(site_index,site_index,t_index);
    end
    sum_diag(t_index) = mean(temp);
    sum_diag_svar(t_index) = mean(temp_svar)/sqrt(NumOfVertexs);
end
errorbar((1:1:TempSlice)*D_Tau,sum_diag,sum_diag_svar);
