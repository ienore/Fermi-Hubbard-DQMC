
%% Set up Parameters
zjy_index = 1;

T_hop = 1.0;
Uene = 4;
%Uene = 0.0;
Beta = 4;
D_Tau = 0.05;

NumOfWarm = 100;
NumOfEpoch = 0;

NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

global N_wrap;
N_wrap = 5;
%% Related Parametres
Miu = Uene/2;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_cut = 5.0;
id_mat = eye(NumOfVertexs);
global expmK;
expmK = expm(D_Tau*(K*T_hop));
global expmK_inv;
expmK_inv = expm(-D_Tau*(K*T_hop));
global B_up_list;
global B_down_list;

%% udv lists
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
global id_mat;

refine_number = fix(TempSlice/N_wrap);
refine_list = zeros([refine_number,1]);
for i = 1:1:refine_number
    refine_list(i) = N_wrap * i;
end
u_up_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);%_1表示前段
d_up_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
v_up_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);

u_down_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
d_down_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
v_down_list_1 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);

u_up_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);%_2表示后端
d_up_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
v_up_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);

u_down_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
d_down_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);
v_down_list_2 = zeros([NumOfVertexs,NumOfVertexs,refine_number]);

%B_up_list = zeros([NumOfVertexs,NumOfVertexs,TempSlice]);
%B_down_list = zeros([NumOfVertexs,NumOfVertexs,TempSlice]);

%% make everything to be identity
for index = 1:1:refine_number
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
for index = 1:1:TempSlice
   B_up_list(:,:,index) = id_mat;
   B_down_list(:,:,index) = id_mat;
end

%% set up measure settings
mea_x = 1;
mea_y = 1;
mea_result_auxi = zeros([1,(TempSlice-1)*NumOfEpoch*2]);
count_list = zeros([1,TempSlice])+1;
count = 1.0;
%% Calculation
for warm_index = 1:1:NumOfWarm
    if mod(zjy_index,8) == 1
        fprintf("D_Tau = %f, Warming_Ratio = %f\n",D_Tau,warm_index/NumOfWarm);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index,N_wrap) == 0 %这个循环是从小更新到大的，所以也是更新list1
                refine_index = time_index/N_wrap;
                Binned_up = id_mat;
                Binned_down = id_mat;
                if time_index >0
                    for B_index = ((refine_index - 1)*N_wrap + 1):1:refine_index*N_wrap
                        Binned_up = B_up_list(:,:,B_index) * Binned_up;
                        Binned_down = B_down_list(:,:,B_index) * Binned_up;
                    end
                end
                [u_up,d_up,v_up] = svdsim(Binned_up);
                [u_down,d_down,v_down] = svdsim(Binned_down);
                u_up_list_1(:,:,refine_index) = u_up;
                d_up_list_1(:,:,refine_index) = d_up;
                v_up_list_1(:,:,refine_index) = v_up;
                u_down_list_1(:,:,refine_index) = u_down;
                d_down_list_1(:,:,refine_index) = d_down;
                v_down_list_1(:,:,refine_index) = v_down;
                green_L_up = Get_G_L_from_list(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L_from_list(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = B_up_list(:,:,time_index);
            B_trans_inv_up = Get_B_L_inv(1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = B_down_list(:,:,time_index);
            B_trans_inv_down = Get_B_L_inv(-1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
            end
        else
            time_index = TempSlice - time_index_mother + 1;
            if  mod(time_index,N_wrap) == 0%这个循环是从大更新到小的，所以更新list2
                refine_index = time_index/N_wrap;
                Binned_up = id_mat;
                Binned_down = id_mat;
                if time_index < TempSlice - N_wrap
                    for B_index = (refine_index *N_wrap + 1):1:(refine_index+1)*N_wrap
                        Binned_up = B_up_list(:,:,B_index) * Binned_up;
                        Binned_down = B_down_list(:,:,B_index) * Binned_up;
                    end
                else
                   for B_index = (time_index + 1):1:TempSlice 
                       Binned_up = B_up_list(:,:,B_index) * Binned_up;
                       Binned_down = B_down_list(:,:,B_index) * Binned_up;
                   end
                end
                [u_up,d_up,v_up] = svdsim(Binned_up);
                [u_down,d_down,v_down] = svdsim(Binned_down);
                u_up_list_2(:,:,refine_index) = u_up;
                d_up_list_2(:,:,refine_index) = d_up;
                v_up_list_2(:,:,refine_index) = v_up;
                u_down_list_2(:,:,refine_index) = u_down;
                d_down_list_2(:,:,refine_index) = d_down;
                v_down_list_2(:,:,refine_index) = v_down;
                green_L_up = Get_G_L_from_list(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L_from_list(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = B_up_list(:,:,time_index+1);
            B_trans_inv_up = Get_B_L_inv(1,time_index+1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = B_down_list(:,:,time_index+1);
            B_trans_inv_down = Get_B_L_inv(-1,time_index+1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
            end
        end %green_L_up/down is calculated done.
        site_mea = randi([1,NumOfVertexs]);
        mea_result_auxi(count) = green_L_up(site_mea,site_mea) * green_L_down(site_mea,site_mea);
        count = count + 1;
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
        B_up_list(:,:,time_index) = Get_B_L(1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        B_down_list(:,:,time_index) = Get_B_L(-1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
    end
    end
end
disp(mean(mea_result_auxi))


