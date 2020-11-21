function  WarmUp(zjy_index,NumOfWarm)
    %һЩ��ע��
    %���������Ȼ�����Ҫ����ʽ�����ؿ��巽������ɹ۲���֮ǰ����
    %Ŀ��������ϵ�ﵽ�趨���¶ȣ����Գ�֮Ϊ�Ȼ�
    %һ�㽨���Ȼ�100�����ϣ��ر�����Ϊlistһ��ʼ�õ�λ�����ʼ���ģ�ǰ����epoch�Ľ����ǳ��ҡ�
    %������Ҫ�Ȼ����٣�����ͨ������һ�¹������ȣ������򲻰����ⲿ�֣����û������Լ���ӡ�
    
    %���������zjy_index ���������Ϊ0��ʱ���Ȼ�����ʾ�Ȼ����ȡ�����Ϊ���м���׼���ġ�
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
    
    green_L_up = id_mat;%�õ�λ��������ʼ��Green����
    green_L_down = id_mat;
    for epoch_index = 1:1: NumOfWarm
        if mod(zjy_index,8) == 0 && mod(epoch_index,NumOfWarm/1000)==0
            fprintf("D_Tau = %f,Warm_Ratio = %f\n",D_Tau,epoch_index/NumOfWarm);
        end
        for reverse_sign = 0:1:1 
            %reverse_sign = 0   time_index = 2,3,4,5,...,TempSlice
            %reverse_sign = 1   time_index = TempSlice-1,TempSlice-2,...,1
        for time_index_mother = 2:1:TempSlice
            %% ���Ƹ���
                %��  G(t+dt,t+dt) = B(t+dt,t)G(t,t)B^(-1)(t+dt,t)  �����Ƶظ���Green����
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
            %% ��ʱ��ľ�ȷ����
                %��һ�δ�����ǰ�ˣ�����ʱ���Ӵ�С��ʱ����udv���������ȷ����Green�õġ�
                %˳ʱ��ľ�ȷ�����ں��棬֮���Էŵ�λ�ò�ͬ������Ϊ��G(t,t)����ʱ���Ƿ��漰t����¡���һϸ�������˳���йص��µ�
                 
            if mod(time_index,N_wrap) == 0 && reverse_sign == 1 %1��ʾ�Ӵ�С��һ��ϸ���ϵĴ�����ǣ�Ҫ��update L��֮ǰ����green_L
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
                [u_up_new,d_up_new,v_up_new] = svdsim(D_up*(V_up'*u_up)*d_up);%��ֵ�ȶ���svd�ֽ�ʱU*D*V'����ԭ���󣬼���V��Ҫת��
                [u_down_new,d_down_new,v_down_new] = svdsim(D_down*(V_down'*u_down)*d_down);
                u_up_list_2(:,:,refine_index) = U_up * u_up_new;
                d_up_list_2(:,:,refine_index) = d_up_new;
                v_up_list_2(:,:,refine_index) = v_up * v_up_new;
                u_down_list_2(:,:,refine_index) = U_down * u_down_new;
                d_down_list_2(:,:,refine_index) = d_down_new;
                v_down_list_2(:,:,refine_index) = v_down * v_down_new;
                green_L_up = Get_G_L_from_list(1.0,time_index);%��ȷ����Green����
                green_L_down = Get_G_L_from_list(-1.0,time_index); 
            end
            %% Sigma�е�Spin��ת
                %���¸�������ÿ��site������һ���Ƿ���£��������Լ�С��������
            for site_index = 1:1:NumOfVertexs
                delta_up = zeros(NumOfVertexs);
                delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
                delta_down = zeros(NumOfVertexs);
                delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
                R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
                R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
                PosToChange = R_up*R_down;%spin��ת�Ľ��ܸ���
                if rand < PosToChange
                    Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                    green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                    green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
                end
            end
            %% ����B����
            B_up_list(:,:,time_index) = Get_B_L(1,time_index);
            B_down_list(:,:,time_index) = Get_B_L(-1,time_index);
            %% ˳ʱ��ľ�ȷ����
                % ÿ��N_wrap�ͼ���һ��B�ۻ���UDV�ֽⲢ����
                %����list�д���1��2�� 1��ʾB(0,t);2��ʾB(Beta,t)
            if mod(time_index,N_wrap) == 0 && reverse_sign == 0 %0��ʾ��С����
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
                [u_up_new,d_up_new,v_up_new] = svdsim(d_up*(v_up'*U_up)*D_up);%��ֵ�ȶ���svd�ֽ�ʱU*D*V'����ԭ���󣬼���V��Ҫת��
                [u_down_new,d_down_new,v_down_new] = svdsim(d_down*(v_down'*U_down)*D_down);
                u_up_list_1(:,:,refine_index) = u_up * u_up_new;
                d_up_list_1(:,:,refine_index) = d_up_new;
                v_up_list_1(:,:,refine_index) = V_up * v_up_new;
                u_down_list_1(:,:,refine_index) = u_down * u_down_new;
                d_down_list_1(:,:,refine_index) = d_down_new;
                v_down_list_1(:,:,refine_index) = V_down * v_down_new;
                if time_index == TempSlice %%Ϊ��ʱ������׼��
                    u_up_list_2(:,:,refine_index) = id_mat;
                    d_up_list_2(:,:,refine_index) = id_mat;
                    v_up_list_2(:,:,refine_index) = id_mat;
                    u_down_list_2(:,:,refine_index) = id_mat;
                    d_down_list_2(:,:,refine_index) = id_mat;
                    v_down_list_2(:,:,refine_index) = id_mat;
                end
                green_L_up = Get_G_L_from_list(1.0,time_index);%��ȷ����Green����
                green_L_down = Get_G_L_from_list(-1.0,time_index);

            end
        end
        end
    end
end