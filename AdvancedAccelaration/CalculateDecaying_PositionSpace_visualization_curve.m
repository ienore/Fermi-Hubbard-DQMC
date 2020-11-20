
%fprintf("FINAL_RUN:%d\n",FINAL_RUN);
zjy_index = 1;
T_hop = 1.0;
NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

Uene = 2;
Miu = Uene/2;
Beta = 2;
D_Tau = 0.2;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
NumOfWarm = 100;
%NumOfWarm = 10;
NumOfEpoch = 100;
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 10;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);

mea_x = 1;
mea_y = 1;
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,(TempSlice-1)*NumOfEpoch*2]);
mea_result_auxi_new = zeros([1,(TempSlice-1)*NumOfEpoch*2]);
mea_result_propa = zeros([TempSlice,2*(TempSlice-1)*NumOfEpoch*NumOfVertexs]);
count_list = zeros([1,TempSlice])+1;
count = 1.0;
count_new = 1;
Sigma = WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
for epoch_index = 1:1:NumOfEpoch
    if mod(zjy_index,8) == 1 && mod(epoch_index,NumOfEpoch/1000)==0
        fprintf("D_Tau = %f,MC_Ratio = %f\n",D_Tau,epoch_index/NumOfEpoch);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index_mother,N_wrap) == 1 || time_index_mother == 2
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
        for site_index_auxi = 1:1:NumOfVertexs
            mea_result_auxi(count) = green_down_c(1,1)+green_up_c(1,1);
            count = count + 1;
        end
        t_start = 0;
        B_propa = eye(NumOfVertexs);
        for propa_time_index = 1:1:TempSlice
            if mod(propa_time_index,N_wrap) == 1 || propa_time_index == 1
                G_propa = Get_G_L2(1,propa_time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
                B_propa = Get_B_L(1,propa_time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                G_propa = B_propa * G_propa;
            end
            
           %B_propa = Get_B_L(1,propa_time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
           %G_propa = B_propa * G_propa;

           for site_index = 1:1:NumOfVertexs
                mea_result_propa(propa_time_index,count_list(propa_time_index)) = G_propa(site_index,site_index);
                count_list(propa_time_index) = count_list(propa_time_index) + 1;
           end
        end
       
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

errorbar((1:1:TempSlice).*D_Tau/Beta,plot_mean_ori,plot_svar_ori,'r');
title(['L = ',num2str(NumInEdge),'  ','Uene = ',num2str(Uene), '   Beta = ',num2str(Beta),'   N_wrap=',num2str(N_wrap)]);
xlabel('Tau ratio');
ylabel('G_{propa}_Err ratio');




