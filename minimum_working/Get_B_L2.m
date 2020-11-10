function B_L2 = Get_B_L2(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2 = eye(NumOfVertexs);
    [U,S,V] = svdsim(B_L2);
    Bin_Size = 4;
    for Bin_index = 1:1:fix((L2-L1)/Bin_Size)+1
        Binned_B_L = eye(NumOfVertexs);
        L_start = L1 + (Bin_index-1) * Bin_Size ;
        if Bin_index < fix((L2-L1)/Bin_Size)+1       
            for index = 1:1:Bin_Size
                Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
            end
        else
            for index = 1:1:(L2 - L_start)
                Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
            end
        end
        [U_new,S_new,V_new] = svdsim(Binned_B_L);
        [u,s,v] = svdsim(S_new*(V_new'*U)*S);
        U = U_new*u;
        S = s;
        V = V*v;
       % [U,S,V_new] = svdsim(Binned_B_L*U*S);
       % V = V * V_new;
    end
    B_L2 = U * S * V';
end

