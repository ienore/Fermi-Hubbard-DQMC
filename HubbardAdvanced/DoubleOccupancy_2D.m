%% 添加子文件夹
addpath(genpath(pwd));
%% Set up Parameters
global NumInEdge;
global NumOfVertexs;
global D_Tau;
global TempSlice;
global Uene;
global Miu;
global Beta;
global lambda;
global T_hop;
global N_wrap;
zjy_index = 0;
T_hop = 1.0;
NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K_2d(NumInEdge);
Uene = 1;
Beta = 4;
D_Tau = 0.05;
NumOfWarm = 100;
NumOfEpoch = 200;
N_wrap = 10;%暂时地，我们要求TempSlice是N_wrap的整数倍
%% Related Parametres
global expmK;
global expmK_inv;
global id_mat;
global Sigma;

TempSlice = Beta/D_Tau;
Miu = Uene/2;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
id_mat = eye(NumOfVertexs);
expmK = expm(D_Tau*(K*T_hop));
expmK_inv = expm(-D_Tau*(K*T_hop));
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
global B_up_list;
global B_down_list;
global u_up_list_1;
global d_up_list_1;
global v_up_list_1;
global u_down_list_1;
global d_down_list_1;
global v_down_list_1;
global u_up_list_2;
global d_up_list_2;
global v_up_list_2;
global u_down_list_2;
global d_down_list_2;
global v_down_list_2;

%% Set up id_mats as initialization
green_L_up = id_mat;
green_L_down = id_mat;
for index = 1:1:TempSlice
   B_up_list(:,:,index) = id_mat;
   B_down_list(:,:,index) = id_mat;
end
for index = 1:1:TempSlice/N_wrap
   u_up_list_1(:,:,index) = id_mat;
   d_up_list_1(:,:,index) = id_mat;
   v_up_list_1(:,:,index) = id_mat;
   u_down_list_1(:,:,index) = id_mat;
   d_down_list_1(:,:,index) = id_mat;
   v_down_list_1(:,:,index) = id_mat;
   u_up_list_2(:,:,index) = id_mat;
   d_up_list_2(:,:,index) = id_mat;
   v_up_list_2(:,:,index) = id_mat;
   u_down_list_2(:,:,index) = id_mat;
   d_down_list_2(:,:,index) = id_mat;
   v_down_list_2(:,:,index) = id_mat;
end

%% Set up Measure Matrixs
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
%%热化
WarmUp(zjy_index,NumOfWarm);
%% Calculation
    %这里大部分程序与WarmUp.m里面一致，请参看里面的注释。
    %为了简洁，我们略去重复的注释。
for epoch_index = 1:1:NumOfEpoch
    %% MonteCarlo前半部分
    if mod(zjy_index,8) == 0 && mod(epoch_index,NumOfEpoch/1000)==0
        fprintf("D_Tau = %f,MC_Ratio = %f\n",D_Tau,epoch_index/NumOfEpoch);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            B_trans_up = B_up_list(:,:,time_index);
            B_trans_inv_up = Get_B_L_inv(1,time_index);
            B_trans_down = B_down_list(:,:,time_index);
            B_trans_inv_down = Get_B_L_inv(-1,time_index);
            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
        else
            time_index = TempSlice - time_index_mother + 1;
            B_trans_up = B_up_list(:,:,time_index+1);
            B_trans_inv_up = Get_B_L_inv(1,time_index+1);
            B_trans_down = B_down_list(:,:,time_index+1);
            B_trans_inv_down = Get_B_L_inv(-1,time_index+1);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
        end
        green_up = green_L_up;
        green_down = green_L_down;
        green_up_c = id_mat - transpose(green_L_up);
        green_down_c = id_mat - transpose(green_L_down);
        energy = 0.0;
        %% 测量可观测量
        for site_index_auxi = 1:1:NumOfVertexs
            mea_result_auxi(count) = green_down_c(site_index_auxi,site_index_auxi)*green_up_c(site_index_auxi,site_index_auxi);
            mea_result_DD(count) = green_down_c(site_index_auxi,site_index_auxi)*green_up_c(site_index_auxi,site_index_auxi);
            energy = energy + Uene * mea_result_DD(count) - Miu * mea_result_auxi(count);
            count = count + 1;
        end
        %% MonteCarlo后半部分
        if mod(time_index,N_wrap) == 0 && reverse_sign == 1 
            refine_index = time_index/N_wrap;
            Binned_up = id_mat;
            Binned_down = id_mat;
            for udv_index = (time_index+1):1:(time_index + N_wrap)
                Binned_up = B_up_list(:,:,udv_index) * Binned_up;
                Binned_down = B_down_list(:,:,udv_index) * Binned_down;
            end
            [u_up,d_up,v_up] = svdsim(Binned_up);
            [u_down,d_down,v_down] = svdsim(Binned_down);
            U_up = u_up_list_2(:,:,refine_index + 1);
            D_up = d_up_list_2(:,:,refine_index + 1);
            V_up = v_up_list_2(:,:,refine_index + 1);
            U_down = u_down_list_2(:,:,refine_index + 1);
            D_down = d_down_list_2(:,:,refine_index + 1);
            V_down = v_down_list_2(:,:,refine_index + 1);
            [u_up_new,d_up_new,v_up_new] = svdsim(D_up*(V_up'*u_up)*d_up);
            [u_down_new,d_down_new,v_down_new] = svdsim(D_down*(V_down'*u_down)*d_down);
            u_up_list_2(:,:,refine_index) = U_up * u_up_new;
            d_up_list_2(:,:,refine_index) = d_up_new;
            v_up_list_2(:,:,refine_index) = v_up * v_up_new;
            u_down_list_2(:,:,refine_index) = U_down * u_down_new;
            d_down_list_2(:,:,refine_index) = d_down_new;
            v_down_list_2(:,:,refine_index) = v_down * v_down_new;
            green_L_up = Get_G_L_from_list(1.0,time_index);
            green_L_down = Get_G_L_from_list(-1.0,time_index);
        end
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
        B_up_list(:,:,time_index) = Get_B_L(1,time_index);
        B_down_list(:,:,time_index) = Get_B_L(-1,time_index);
        if mod(time_index,N_wrap) == 0 && reverse_sign == 0
            refine_index = time_index/N_wrap;
            Binned_up = id_mat;
            Binned_down = id_mat;
            for udv_index = (time_index-N_wrap+1):1:time_index
                Binned_up = B_up_list(:,:,udv_index) * Binned_up;
                Binned_down = B_down_list(:,:,udv_index) * Binned_down;
            end
            [u_up,d_up,v_up] = svdsim(Binned_up);
            [u_down,d_down,v_down] = svdsim(Binned_down);
            if time_index > N_wrap
               U_up = u_up_list_1(:,:,refine_index - 1);
               D_up = d_up_list_1(:,:,refine_index - 1);
               V_up = v_up_list_1(:,:,refine_index - 1);
               U_down = u_down_list_1(:,:,refine_index - 1);
               D_down = d_down_list_1(:,:,refine_index - 1);
               V_down = v_down_list_1(:,:,refine_index - 1);
            else
                U_up = id_mat;
                D_up = id_mat;
                V_up = id_mat;
                U_down = id_mat;
                D_down = id_mat;
                V_down = id_mat;
            end
            [u_up_new,d_up_new,v_up_new] = svdsim(d_up*(v_up'*U_up)*D_up);
            [u_down_new,d_down_new,v_down_new] = svdsim(d_down*(v_down'*U_down)*D_down);
            u_up_list_1(:,:,refine_index) = u_up * u_up_new;
            d_up_list_1(:,:,refine_index) = d_up_new;
            v_up_list_1(:,:,refine_index) = V_up * v_up_new;
            u_down_list_1(:,:,refine_index) = u_down * u_down_new;
            d_down_list_1(:,:,refine_index) = d_down_new;
            v_down_list_1(:,:,refine_index) = V_down * v_down_new;
            if time_index == TempSlice 
                u_up_list_2(:,:,refine_index) = id_mat;
                d_up_list_2(:,:,refine_index) = id_mat;
                v_up_list_2(:,:,refine_index) = id_mat;
                u_down_list_2(:,:,refine_index) = id_mat;
                d_down_list_2(:,:,refine_index) = id_mat;
                v_down_list_2(:,:,refine_index) = id_mat;
            end
            green_L_up = Get_G_L_from_list(1.0,time_index);
            green_L_down = Get_G_L_from_list(-1.0,time_index);
        end   
    end
    end
end


%% 最后数据的处理
    %有两个函数值得推荐：
    %Bin.m 返回平均值与误差
    %Calc_Correlation.m 返回关联函数
sample = mea_result_auxi;
[mean_sample,svar_sample] = Bin(sample);
[step_list,ans_list] = Calc_Correlation(sample);
fprintf("Mean=%f\t svar=%f\t error_ratio=%f\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));
plot(step_list,ans_list);


