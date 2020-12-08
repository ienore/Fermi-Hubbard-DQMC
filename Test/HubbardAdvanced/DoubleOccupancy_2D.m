%% 核心公式：Green函数
% G_ij(t,t) = <C_i(t) C_j^\dagger(t)>
% G_ij(t,t) = inv(I + B(t,0)B(\beta,t))

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
NumOfVertexs = NumInEdge^2;%维度，如果2维的话，SNN的表达式也许要改变
K = Get_K_2d(NumInEdge);%维度
Uene = 2.0;
Beta = 2;
D_Tau = 0.05;
NumOfWarm = 50;
NumOfEpoch = 50;
N_wrap = 10;%暂时地，我们要求TempSlice是N_wrap的整数倍
%% 与任务相关的新定义
mea_gap = 5;%采样间隔
num_mea = TempSlice * 2 * NumOfEpoch / mea_gap;%采样上限
S_px = pi;
S_py = pi;
mea_px = pi/2;
mea_py = pi/2;
eta = 0.037;
mea_result_SNN = zeros([1,num_mea]);
mea_result_E = zeros([1,num_mea]);
mea_result_DD = zeros([1,num_mea]);
mea_result_S = zeros([1,num_mea]);
mea_result_green_mat = zeros([NumOfVertexs,NumOfVertexs,num_mea]);
mea_result_green_tau_mat = zeros([NumOfVertexs,NumOfVertexs,TempSlice,num_mea]);
mea_result_Gb2_mat = zeros([NumOfVertexs,NumOfVertexs,num_mea]);
mea_result_gb2 = zeros([1,num_mea]);
count_mea = 0;

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
%辅助测量矩阵
mea_result_auxi = [];
%% 热化
[green_L_up,green_L_down] = WarmUp(zjy_index,NumOfWarm);
%% Calculation
    %这里大部分程序与WarmUp.m里面一致，请参看里面的注释。
    %为了简洁，我们略去重复的注释。
for epoch_index = 1:1:NumOfEpoch
    time_choice = randi([1,TempSlice]);%随机选一个时刻做测量
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
        energy_k = 0.0;
        
        %% 测量可观测量
        if (mod(abs(time_index-time_choice),mea_gap)==0 )
            count_mea = count_mea + 1;%测量指标
            %% 双占据
            for site_index_auxi = 1:1:NumOfVertexs
                mea_result_auxi(count_mea) = green_down_c(site_index_auxi,site_index_auxi)*green_up_c(site_index_auxi,site_index_auxi);
            end
            mea_result_DD(count_mea) = mean(diag(green_down_c).*diag(green_up_c));
            %% 能量
            for site_index = 1:1:NumOfVertexs
                energy_k = energy_k + Uene * (green_down_c(site_index,site_index))*(green_up_c(site_index,site_index)) ;
            end
            for index_1 = 1:1:NumOfVertexs
                for index_2 = 1:1:NumOfVertexs
                    if K(index_1,index_2) ~= 0
                        energy_k = energy_k - T_hop *(green_down_c(index_1,index_2) + green_up_c(index_1,index_2));
                    end
                end
            end
            mea_result_E(count_mea) = energy_k/(NumOfVertexs);
           %% 最近邻S
            SNN = 0.0;
            a = 0.0;
            for site_i = 1:1:NumOfVertexs
                for site_j = 1:1:NumOfVertexs
                    if K(site_i,site_j) ~= 0
                        a = green_up_c(site_i,site_i) * green_up_c(site_j,site_j) + green_up_c(site_i,site_j) * green_up(site_i,site_j) + ...
                        green_down_c(site_i,site_i) * green_down_c(site_j,site_j) + green_down_c(site_i,site_j) * green_down(site_i,site_j) - ...
                        green_down_c(site_i,site_i) * green_up_c(site_j,site_j) - green_up_c(site_i,site_i) * green_down_c(site_j,site_j);
                        SNN = SNN + a/(4*4);%4来自自旋，6来自最近邻，如果2d就要变成4
                    end
                end
            end
            SNN = real(SNN)/NumOfVertexs;
            mea_result_SNN(count_mea) = SNN ;
           %% 结构因子
            S = 0.0;
            S_px = pi;
            S_py = pi;
            for site_i = 1:1:NumOfVertexs
                for site_j = 1:1:NumOfVertexs
                    [xi,yi] = IndexToCoor_2d(site_i,NumInEdge);
                    [xj,yj] = IndexToCoor_2d(site_j,NumInEdge);
                    factor = exp(1i*(S_px*(xi-xj)+S_py*(yi-yj)));

                    a = green_up_c(site_i,site_i) * green_up_c(site_j,site_j) + green_up_c(site_i,site_j) * green_up(site_i,site_j) + ...
                    green_down_c(site_i,site_i) * green_down_c(site_j,site_j) + green_down_c(site_i,site_j) * green_down(site_i,site_j) - ...
                    green_down_c(site_i,site_i) * green_up_c(site_j,site_j) - green_up_c(site_i,site_i) * green_down_c(site_j,site_j);

                    S = S + factor * a/(NumOfVertexs);
                end
            end
            mea_result_S(count_mea) = real(S) ;
           %% 测量等时Green函数
%             mea_result_green_mat(:,:,count_mea) = green_L_up + green_L_down;
           %% 测量非等时Green函数
%            for mea_time_index = 1:1:TempSlice
%                 if mod(mea_time_index,N_wrap) == 0 || mea_time_index <= N_wrap
%                     mea_G_L2_up = Get_G_L2_old(1,mea_time_index);
%                     mea_G_L2_down = Get_G_L2_old(-1,mea_time_index);
%                 else
%                     mea_G_L2_up = B_up_list(:,:,mea_time_index) * mea_G_L2_up;
%                     mea_G_L2_down = B_down_list(:,:,mea_time_index) * mea_G_L2_down;
%                 end
%                 mea_result_green_tau_mat(:,:,mea_time_index,count_mea) = mea_G_L2_up + mea_G_L2_down;
%            end
           %% 测量G(beta/2)
            mea_result_Gb2_mat(:,:,count_mea) = Get_G_L2_old(1,TempSlice/2)+Get_G_L2_old(2,TempSlice/2);
            mea_result_gb2(count_mea) = mean(diag(mea_result_Gb2_mat(:,:,count_mea)));
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

% %% 处理Green函数
% %等时部分
% green_mat = zeros(NumOfVertexs);
% for index = 1:1:count_mea
%     green_mat = green_mat + mea_result_green_mat(:,:,index)/count_mea;
% end
% %半时部分
% 
% 
% %非等时部分
% green_mat_list = zeros([NumOfVertexs,NumOfVertexs,TempSlice]);
% for time_index = 1:1:TempSlice
%    for mea_index = 1:1:count_mea
%       green_mat_list(:,:,time_index) = green_mat_list(:,:,time_index) + ...
%           mea_result_green_tau_mat(:,:,time_index,mea_index)/count_mea;
%    end
% end
% %% 特别处理（与当前任务相关）
% mea_px = 0;
% mea_py = 0;
% G_k_tau = zeros([1,TempSlice]);
% for time_index = 1:1:TempSlice
%     temp = 0.0;
%     for site_i = 1:1:NumOfVertexs
%         for site_j = 1:1:NumOfVertexs
%             [xi,yi] = IndexToCoor_2d(site_i,NumInEdge);
%             [xj,yj] = IndexToCoor_2d(site_j,NumInEdge);
%             factor = exp(1i*(mea_px*(xi-xj)+mea_py*(yi-yj)));
%             temp = temp + factor * green_mat_list(site_i,site_j,time_index);
%         end
%     end
%     temp = real(temp)/NumOfVertexs;
%     G_k_tau(time_index) = temp;
% end
% plot((1:1:TempSlice)*D_Tau,G_k_tau);
%% 最后数据的处理
    %有两个函数值得推荐：
    %Bin.m 返回平均值与误差
    %Calc_Correlation.m 返回关联函数
    
sample = mea_result_E(1:count_mea);
[mean_sample,svar_sample] = Bin(sample);
%[step_list,ans_list] = Calc_Correlation(sample);
fprintf("E_k : Mean=%f\t svar=%e\t error_ratio=%e\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));

sample = mea_result_DD(1:count_mea);
[mean_sample,svar_sample] = Bin(sample);
%[step_list,ans_list] = Calc_Correlation(sample);
fprintf("DD : Mean=%f\t svar=%e\t error_ratio=%e\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));

sample = mea_result_S(1:count_mea);
[mean_sample,svar_sample] = Bin(sample);
%[step_list,ans_list] = Calc_Correlation(sample);
fprintf("S : Mean=%f\t svar=%e\t error_ratio=%e\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));

sample = mea_result_SNN(1:count_mea);
[mean_sample,svar_sample] = Bin(sample);
%[step_list,ans_list] = Calc_Correlation(sample);
fprintf("SNN : Mean=%f\t svar=%e\t error_ratio=%e\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));

sample = mea_result_gb2(1:count_mea)*Beta;
[mean_sample,svar_sample] = Bin(sample);
%[step_list,ans_list] = Calc_Correlation(sample);
fprintf("GB2: Mean=%f\t svar=%e\t error_ratio=%e\n",mean_sample,svar_sample,abs(svar_sample/mean_sample));
