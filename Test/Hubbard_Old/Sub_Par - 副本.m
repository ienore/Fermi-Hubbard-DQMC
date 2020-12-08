function  [return_mean,return_svar] = Sub_Par(FINAL_RUN,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop)
%fprintf("FINAL_RUN:%d\n",FINAL_RUN);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WarmUp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
    for epoch = 1:1:NumOfWarm
        %disp(epoch);
        O_Plus = Get_O_alpha(1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        O_Minus = Get_O_alpha(-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        Green_Plus_T = inv(O_Plus);
        Green_Minus_T = inv(O_Minus);
        for l_change = 1:1:TempSlice
            for i_change = 1:1:NumOfVertexs


                g_Plus = Green_Plus_T((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change);
                g_Minus = Green_Minus_T((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change);
                p_Plus = 1.0+(1.0-g_Plus)*(exp(-2*lambda*Sigma(l_change,i_change))-1.0);
                p_Minus = 1.0+(1.0-g_Minus)*(exp(2*lambda*Sigma(l_change,i_change))-1.0);
                p = p_Minus*p_Plus;
                PosToChange = p/(1.0+p);

                if ((rand < PosToChange))
                    t_Plus = (exp(-2*lambda*Sigma(l_change,i_change))-1.0)/(p_Plus);
                    t_Plus_Mat = zeros([TempSlice*NumOfVertexs,TempSlice*NumOfVertexs]);
                    t_Plus_Mat((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change) = t_Plus;
                    Green_Plus_T = Green_Plus_T + (Green_Plus_T - eye(TempSlice*NumOfVertexs))*t_Plus_Mat*Green_Plus_T;

                    t_Minus = (exp(2*lambda*Sigma(l_change,i_change))-1.0)/(p_Minus);
                    t_Minus_Mat = zeros([TempSlice*NumOfVertexs,TempSlice*NumOfVertexs]);
                    t_Minus_Mat((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change) = t_Minus;
                    Green_Minus_T = Green_Minus_T + (Green_Minus_T - eye(TempSlice*NumOfVertexs))*t_Minus_Mat*Green_Minus_T;
                    Sigma(l_change,i_change) = - Sigma(l_change,i_change);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%PhysicalCalculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    fprintf("Warmed\n");
    big_temp_1 = zeros([NumOfEpoch*TempSlice,1]);%%%%%%%%%%%%%%%%%%%%%555555
    big_temp_2 = zeros([NumOfEpoch*NumOfVertexs*TempSlice,1]);
    temp_observable = zeros([NumOfEpoch*NumOfVertexs*TempSlice,1]);
    count = 1;
    for i = 1:1:NumOfEpoch
        l_change = randi(TempSlice);
        i_change = randi(NumOfVertexs);

        if( mod(i,10000)<10 )
            O_Plus = Get_O_alpha(1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            O_Minus = Get_O_alpha(-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            Green_Plus_T = inv(O_Plus);
            Green_Minus_T = inv(O_Minus);
        end
        green_L_up = Get_G_L(1.0,l_change,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        green_L_down = Get_G_L(-1.0,l_change,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        g_Plus = green_L_up(i_change,i_change);
        g_Minus = green_L_down(i_change,i_change);
        p_Plus = 1.0+(1.0-g_Plus)*(exp(-2*lambda*Sigma(l_change,i_change))-1.0);
        p_Minus = 1.0+(1.0-g_Minus)*(exp(2*lambda*Sigma(l_change,i_change))-1.0);
        p = p_Minus*p_Plus;
        PosToChange = p/(1.0+p);

        if ((rand < PosToChange))
            t_Plus = (exp(-2*lambda*Sigma(l_change,i_change))-1.0)/(p_Plus);
            t_Plus_Mat = zeros([TempSlice*NumOfVertexs,TempSlice*NumOfVertexs]);
            t_Plus_Mat((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change) = t_Plus;
            Green_Plus_T = Green_Plus_T + (Green_Plus_T - eye(TempSlice*NumOfVertexs))*t_Plus_Mat*Green_Plus_T;

            t_Minus = (exp(2*lambda*Sigma(l_change,i_change))-1.0)/(p_Minus);
            t_Minus_Mat = zeros([TempSlice*NumOfVertexs,TempSlice*NumOfVertexs]);
            t_Minus_Mat((l_change-1)*NumOfVertexs+i_change,(l_change-1)*NumOfVertexs+i_change) = t_Minus;
            Green_Minus_T = Green_Minus_T + (Green_Minus_T - eye(TempSlice*NumOfVertexs))*t_Minus_Mat*Green_Minus_T;
            Sigma(l_change,i_change) = - Sigma(l_change,i_change);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% measure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for L = 1:1:TempSlice
            green_up = Green_Plus_T((L-1)*NumOfVertexs+1:L*NumOfVertexs,(L-1)*NumOfVertexs+1:L*NumOfVertexs);
            green_down = Green_Minus_T((L-1)*NumOfVertexs+1:L*NumOfVertexs,(L-1)*NumOfVertexs+1:L*NumOfVertexs);
            green_up_c = eye(NumOfVertexs) - transpose(green_up);
            green_down_c = eye(NumOfVertexs) - transpose(green_down);
            sum_measure = 0.0;
            for site_i = 1:1:NumOfVertexs
                for site_j = 1:1:NumOfVertexs
                   factor = 0;
                   [i_x,i_y] = IndexToCoor(site_i,NumInEdge);
                   [j_x,j_y] = IndexToCoor(site_j,NumInEdge);
                   delta_r = abs(i_x-j_x) + abs(i_y-j_y);
                   if mod(delta_r,2) == 0
                       factor = 1.0;
                   else
                       factor = -1.0;
                   end
                   a_measure = green_up_c(site_i,site_i) * green_up_c(site_j,site_j) + green_up_c(site_i,site_j) * green_up(site_i,site_j) + ...
                               green_down_c(site_i,site_i) * green_down_c(site_j,site_j) + green_down_c(site_i,site_j) * green_down(site_i,site_j) - ...
                               green_down_c(site_i,site_i) * green_up_c(site_j,site_j) - green_up_c(site_i,site_i) * green_down_c(site_j,site_j);
                   sum_measure = sum_measure + factor*a_measure;
                end
            end
            sum_measure = sum_measure / (4*NumOfVertexs*NumOfVertexs);
            big_temp_1(count) = sum_measure;
            count = count +1.0;
        end

        %%%%%%%%%%%%%%%%%%%%%%% measure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    temp_observable = big_temp_1;
    return_mean = mean(temp_observable);
    return_svar = sqrt(var(temp_observable));
end

