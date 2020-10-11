
%fprintf("FINAL_RUN:%d\n",FINAL_RUN);
zjy_index = 1;
T_hop = 1.0;
NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

Uene = 8.0;
Miu = Uene/2;
Beta = 8;
D_Tau = 0.2;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
NumOfWarm = 100;
%NumOfWarm = 0;
NumOfEpoch = 100;
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 5;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);

mea_x = 1;
mea_y = 1;
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_propa = zeros([TempSlice,2*(TempSlice-1)*NumOfEpoch*NumOfVertexs]);
count_list = zeros([1,TempSlice])+1;
t_range = 0.0*mea_result_auxi;
count = 1.0;

WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
for epoch_index = 1:1:NumOfEpoch
    if mod(zjy_index,8) == 1 && mod(epoch_index,NumOfEpoch/1000)==0
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE ?? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        green_up = green_L_up;
        green_down = green_L_down;
        green_up_c = id_mat - transpose(green_L_up);
        green_down_c = id_mat - transpose(green_L_down);
        %green_up = Get_G_L(1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        %green_down = Get_G_L(-1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        
        mea_result_auxi(count) = green_down(1,1)+green_up(1,1);
        t_range(count) = time_index/TempSlice;
        count = count + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE ?? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot_mean = zeros([1,TempSlice]);
plot_svar = zeros([1,TempSlice]);
plot_mean_ori = zeros([1,TempSlice]);
plot_svar_ori = zeros([1,TempSlice]);
for time_index = 1:1:TempSlice
    sample = real(log(mea_result_propa(time_index,:)));
    [plot_mean(time_index),plot_svar(time_index)] = Bin(sample);
end
for time_index = 1:1:TempSlice
    sample_ori = real((mea_result_propa(time_index,:)));
    [plot_mean_ori(time_index),plot_svar_ori(time_index)] = Bin(sample_ori);
end

% errorbar((1:1:TempSlice).*D_Tau/Beta,plot_mean,plot_svar,'r');
% title(['Uene = ',num2str(Uene), '   Beta = ',num2str(Beta)]);
% xlabel('Tau ratio');
% ylabel('log(G_propa)');
% 
% errorbar((1:1:TempSlice).*D_Tau/Beta,plot_mean_ori,plot_svar_ori);
% title(['L = ',num2str(NumInEdge),'  ','Uene = ',num2str(Uene), '   Beta = ',num2str(Beta)]);
% xlabel('Tau ratio');
% ylabel('G_{propa}_Err ratio');

plot(mea_result_auxi);
hold on
plot(t_range)

% paper_mean = zeros([1,TempSlice]);
% paper_svar = zeros([1,TempSlice]);
% for time_index = 1:1:TempSlice
%     temp_data = mea_result_propa(1,:) .* mea_result_propa(time_index,:);
%     [paper_mean(time_index),paper_svar(time_index)] = JackKnife(temp_data,length(temp_data));
% end
% errorbar((1:1:TempSlice).*D_Tau,paper_mean,paper_svar,'r');
    