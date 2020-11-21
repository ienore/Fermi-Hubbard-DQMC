function  WarmUp(zjy_index,NumOfWarm)
    %一些备注：
    %本函数是热化，需要在正式用蒙特卡洛方法计算可观测量之前调用
    %目的是让体系达到设定的温度，所以称之为热化
    %一般建议热化100次以上，特别是因为list一开始用单位矩阵初始化的，前几个epoch的结果会非常乱。
    %具体需要热化多少，可以通过计算一下关联长度，本程序不包含这部分，但用户可以自己添加。
    
    %输入参数：zjy_index 当这个参数为0的时候，热化会显示热化进度。这是为并行计算准备的。
    global NumOfVertexs;
    global D_Tau;
    global TempSlice;
    global lambda;
    global N_wrap;
    global id_mat;
    global Sigma;
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
    
    green_L_up = id_mat;%用单位矩阵来初始化Green函数
    green_L_down = id_mat;
    for epoch_index = 1:1: NumOfWarm
        if mod(zjy_index,8) == 0 && mod(epoch_index,NumOfWarm/1000)==0
            fprintf("D_Tau = %f,Warm_Ratio = %f\n",D_Tau,epoch_index/NumOfWarm);
        end
        for reverse_sign = 0:1:1 
            %reverse_sign = 0   time_index = 2,3,4,5,...,TempSlice
            %reverse_sign = 1   time_index = TempSlice-1,TempSlice-2,...,1
        for time_index_mother = 2:1:TempSlice
            %% 近似更新
                %用  G(t+dt,t+dt) = B(t+dt,t)G(t,t)B^(-1)(t+dt,t)  来近似地更新Green函数
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
            %% 逆时间的精确更新
                %这一段代码提前了，是逆时（从大到小）时更新udv矩阵进而精确更新Green用的。
                %顺时间的精确更新在后面，之所以放的位置不同，是因为算G(t,t)函数时“是否涉及t层更新”这一细节与更新顺序有关导致的
                 
            if mod(time_index,N_wrap) == 0 && reverse_sign == 1 %1表示从大到小，一个细节上的处理就是，要在update L层之前计算green_L
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
                [u_up_new,d_up_new,v_up_new] = svdsim(D_up*(V_up'*u_up)*d_up);%数值稳定，svd分解时U*D*V'才是原矩阵，即，V需要转置
                [u_down_new,d_down_new,v_down_new] = svdsim(D_down*(V_down'*u_down)*d_down);
                u_up_list_2(:,:,refine_index) = U_up * u_up_new;
                d_up_list_2(:,:,refine_index) = d_up_new;
                v_up_list_2(:,:,refine_index) = v_up * v_up_new;
                u_down_list_2(:,:,refine_index) = U_down * u_down_new;
                d_down_list_2(:,:,refine_index) = d_down_new;
                v_down_list_2(:,:,refine_index) = v_down * v_down_new;
                green_L_up = Get_G_L_from_list(1.0,time_index);%精确更新Green函数
                green_L_down = Get_G_L_from_list(-1.0,time_index); 
            end
            %% Sigma中的Spin翻转
                %更新辅助场，每个site都计算一下是否更新，这样可以减小关联长度
            for site_index = 1:1:NumOfVertexs
                delta_up = zeros(NumOfVertexs);
                delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
                delta_down = zeros(NumOfVertexs);
                delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
                R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
                R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
                PosToChange = R_up*R_down;%spin翻转的接受概率
                if rand < PosToChange
                    Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                    green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                    green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
                end
            end
            %% 更新B矩阵
            B_up_list(:,:,time_index) = Get_B_L(1,time_index);
            B_down_list(:,:,time_index) = Get_B_L(-1,time_index);
            %% 顺时间的精确更新
                % 每隔N_wrap就计算一次B累积的UDV分解并保存
                %其中list有代号1和2： 1表示B(0,t);2表示B(Beta,t)
            if mod(time_index,N_wrap) == 0 && reverse_sign == 0 %0表示从小到大
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
                [u_up_new,d_up_new,v_up_new] = svdsim(d_up*(v_up'*U_up)*D_up);%数值稳定，svd分解时U*D*V'才是原矩阵，即，V需要转置
                [u_down_new,d_down_new,v_down_new] = svdsim(d_down*(v_down'*U_down)*D_down);
                u_up_list_1(:,:,refine_index) = u_up * u_up_new;
                d_up_list_1(:,:,refine_index) = d_up_new;
                v_up_list_1(:,:,refine_index) = V_up * v_up_new;
                u_down_list_1(:,:,refine_index) = u_down * u_down_new;
                d_down_list_1(:,:,refine_index) = d_down_new;
                v_down_list_1(:,:,refine_index) = V_down * v_down_new;
                if time_index == TempSlice %%为逆时更新做准备
                    u_up_list_2(:,:,refine_index) = id_mat;
                    d_up_list_2(:,:,refine_index) = id_mat;
                    v_up_list_2(:,:,refine_index) = id_mat;
                    u_down_list_2(:,:,refine_index) = id_mat;
                    d_down_list_2(:,:,refine_index) = id_mat;
                    v_down_list_2(:,:,refine_index) = id_mat;
                end
                green_L_up = Get_G_L_from_list(1.0,time_index);%精确更新Green函数
                green_L_down = Get_G_L_from_list(-1.0,time_index);

            end
        end
        end
    end
end