function parellel_return = CalculateDecaying_ParellelKernel(zjy_index,Beta,D_Tau,NumOfWarm,NumOfWarm_inside,N_wrap,Uene,T_hop,NumInEdge,OutputNumber,NumOfEpochInBin)

    NumOfVertexs = NumInEdge^2;
    TempSlice = Beta/D_Tau;
    lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
    Miu = Uene/2;
    K = Get_K(NumInEdge);
    Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
    id_mat = eye(NumOfVertexs);

    mea_result_auxi = zeros([1,(TempSlice-1)*NumOfEpochInBin*2]);
    parellel_return = zeros([NumOfVertexs,NumOfVertexs,TempSlice,OutputNumber]);

    WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpochInBin,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
    for output_index = 1:1:OutputNumber
        WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm_inside,NumOfEpochInBin,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
        count_list = zeros([1,TempSlice])+1;
        count = 1.0;
        for epoch_index = 1:1:NumOfEpochInBin
            if mod(zjy_index,8) == 1
                fprintf("D_Tau = %f,MC_Ratio = %f, Bin_index = %d\n",D_Tau,epoch_index/NumOfEpochInBin,output_index);
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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE ¡ý %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                green_up = green_L_up;
                green_down = green_L_down;
                green_up_c = id_mat - transpose(green_L_up);
                green_down_c = id_mat - transpose(green_L_down);
                mea_result_auxi(count) = green_down_c(1,1)+green_up_c(1,1);
                count = count + 1;
                for propa_time_index = 1:1:TempSlice
                    if mod(propa_time_index,N_wrap) == 1
                        G_propa = Get_G_L2(1,propa_time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                    else
                        B_propa = Get_B_L(1,propa_time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                        G_propa = B_propa * G_propa;
                    end
                   parellel_return(:,:,propa_time_index,output_index) = parellel_return(:,:,propa_time_index,output_index) + G_propa/(2*(TempSlice-1)*NumOfEpochInBin);
                   count_list(propa_time_index) = count_list(propa_time_index) + 1;
                end

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
    end

end

