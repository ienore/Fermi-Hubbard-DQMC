%fprintf("FINAL_RUN:%d\n",FINAL_RUN);
zjy_index = 1;
T_hop = 1.0;
NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

Uene = 3.0;
Miu = Uene/2;
Beta = 2;
D_Tau = 0.125;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
NumOfWarm = 50;
NumOfWarm_inside = 0;
%NumOfWarm = 0;
NumOfEpoch = 100;
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 10;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);

mea_px = pi;
mea_py = 0;
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_propa = zeros([TempSlice,(TempSlice-1)*NumOfEpoch*2]);
count = 1.0;

WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
for epoch_index = 1:1:NumOfEpoch
    WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm_inside,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);%reduce the correlation
    if mod(zjy_index,8) == 1
        fprintf("D_Tau = %f,MC_Ratio = %f\n",D_Tau,epoch_index/NumOfEpoch);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index_mother,N_wrap) == 2
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
            end
        else
            time_index = TempSlice - time_index_mother + 1;
            if mod(time_index_mother,N_wrap) == 1
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE ¡ý %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        green_up = green_L_up;
        green_down = green_L_down;
        green_up_c = id_mat - transpose(green_L_up);
        green_down_c = id_mat - transpose(green_L_down);
        mea_result_auxi(count) = green_down_c(1,1)*green_up_c(1,1);
        t_start = 0;
        G_equal = Get_G_L(1.0,t_start,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);%up state
        for propa_time_index = 1:1:TempSlice
           B_propa = Get_B_L2(1,propa_time_index,t_start,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene); 
           G_propa = B_propa * G_equal;
           %G_propa = id_mat - transpose(G_propa);
           temp_sum = 0.0;
           for site_m = 1:1:NumOfVertexs
               for site_n = 1:1:NumOfVertexs
                    [mx,my] = IndexToCoor(site_m,NumInEdge);
                    [nx,ny] = IndexToCoor(site_n,NumInEdge);
                    delta_x = mx - nx;
                    delta_y = my - ny;
                    temp_sum = temp_sum + real(exp(1i * (mea_px*delta_x + mea_py*delta_y))*G_propa(site_m,site_n));
               end
           end
           temp_sum = temp_sum / NumOfVertexs;
           mea_result_propa(propa_time_index,count) = temp_sum;
        end
        count = count + 1;
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE ¡ü %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for site_index = 1:1:NumOfVertexs
            delta_up = zeros(NumOfVertexs);
            delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
            delta_down = zeros(NumOfVertexs);
            delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
            R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
            R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
            PosToChange = R_up*R_down;
            if rand < PosToChange
                Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
            end
        end
    end
    end
end
propa_mean = zeros([1,TempSlice]);
propa_svar = zeros([1,TempSlice]);
for time_index = 1:1:TempSlice
    sample = mea_result_propa(time_index,:);
    [propa_mean(time_index),propa_svar(time_index)] = Bin(sample);
end
errorbar((1:1:TempSlice).*D_Tau/Beta,log(propa_mean),propa_svar/propa_mean*log(propa_mean),'r');
title(['px = ',num2str(mea_px), '  py = ',num2str(mea_py),'   Uene = ',num2str(Uene), '   Beta = ',num2str(Beta)]);
xlabel('Tau ratio');
ylabel('log G_propa');

p = polyfit((1:1:TempSlice).*D_Tau,log(propa_mean),1);
fprintf("The fitted slope = %f\n",p(1))
% paper_mean = zeros([1,TempSlice]);
% paper_svar = zeros([1,TempSlice]);
% for time_index = 1:1:TempSlice
%     temp_data = mea_result_propa(1,:) .* mea_result_propa(time_index,:);
%     [paper_mean(time_index),paper_svar(time_index)] = JackKnife(temp_data,length(temp_data));
% end
% errorbar((1:1:TempSlice).*D_Tau,paper_mean,paper_svar,'r');
    
