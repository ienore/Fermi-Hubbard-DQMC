
%% Set up Parameters
zjy_index = 1;
T_hop = 1.0;
NumInEdge = 4;
NumOfVertexs = NumInEdge^3;
K = Get_K_3d(NumInEdge);

Uene = 4;

Beta = 5;
D_Tau = 0.1;
TempSlice = Beta/D_Tau;
NumOfWarm = 1;
%NumOfWarm = 10;
NumOfEpoch = 10;

%% Related Parametres
Miu = Uene/2;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 10;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);

mea_x = 1;
mea_y = 1;
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,(TempSlice-1)*NumOfEpoch*NumOfVertexs*2]);
mea_result_DD = zeros([1,(TempSlice-1)*NumOfEpoch*NumOfVertexs*2]);
mea_result_E = zeros([1,(TempSlice-1)*NumOfEpoch*2]);
count_list = zeros([1,TempSlice])+1;
count = 1.0;
count_new = 1;
count_E = 1;
%% Calculation
WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
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
        energy = 0.0;
        for site_index_auxi = 1:1:NumOfVertexs
            mea_result_auxi(count) = green_down_c(site_index_auxi,site_index_auxi)+green_up_c(site_index_auxi,site_index_auxi);
            mea_result_DD(count) = green_down_c(site_index_auxi,site_index_auxi)*green_up_c(site_index_auxi,site_index_auxi);
            energy = energy + Uene * mea_result_DD(count) - Miu * mea_result_auxi(count);
            count = count + 1;
        end
        energy = 0.0;
        for index_1 = 1:1:NumOfVertexs
            for index_2 = 1:1:NumOfVertexs
                if K(index_1,index_2) ~= 0
                    energy = energy - T_hop * (green_down_c(index_1,index_2) + green_up_c(index_1,index_2));
                end
            end
        end
        mea_result_E(count_E) = energy/(NumOfVertexs);
        count_E = count_E + 1;
       
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


%% Plot the result 
%plot(mea_result_auxi)
[mean_DD,svar_DD] = Bin(mea_result_DD);
%errorbar(mean_DD,svar_DD)
%hold on
%errorbar(0.1554,0,'r*')
disp(mean(mea_result_E/2))



