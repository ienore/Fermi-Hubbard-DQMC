function  [return_mean,return_svar] = Sub_Par2(FINAL_RUN,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop)
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
    big_temp_1 = zeros([NumOfEpoch*TempSlice*NumOfVertexs,1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
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
        
        %%%%%%%%%%%%%%%%%%%%%%% measure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for L = 1:1:TempSlice
            green_up = Green_Plus_T((L-1)*NumOfVertexs+1:L*NumOfVertexs,(L-1)*NumOfVertexs+1:L*NumOfVertexs);
            green_down = Green_Minus_T((L-1)*NumOfVertexs+1:L*NumOfVertexs,(L-1)*NumOfVertexs+1:L*NumOfVertexs);
            green_up_c = eye(NumOfVertexs) - transpose(green_up);
            green_down_c = eye(NumOfVertexs) - transpose(green_down);
            for site_i = 1:1:NumOfVertexs
                big_temp_1(count) = green_up_c(site_i,site_i)*green_down_c(site_i,site_i);
                count = count +1.0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%% measure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    temp_observable = big_temp_1;
    return_mean = mean(temp_observable);
    return_svar = sqrt(var(temp_observable));
end

