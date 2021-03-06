function [U,S,V] = Get_B_L2_inv_svd(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2_inv = eye(NumOfVertexs);
    [U,S,V] = svdsim(B_L2_inv);
    Bin_Size = 7;
    
    if fix((L2-L1)/Bin_Size)==(L2-L1)/Bin_Size
        Bin_run = (L2-L1)/Bin_Size;
    else
        Bin_run = fix((L2-L1)/Bin_Size)+1;
    end
    
    for Bin_index = 1:1:Bin_run
        Binned_B_L = eye(NumOfVertexs);
        L_start = L1 + (Bin_index-1) * Bin_Size ;
        if Bin_index < Bin_run    
            for index = 1:1:Bin_Size
                Binned_B_L = Binned_B_L*Get_B_L_inv(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            end
        else
            for index = 1:1:(L2 - L_start)
                Binned_B_L = Binned_B_L*Get_B_L_inv(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            end
        end
        [U_new,S_new,V_new] = svdsim(Binned_B_L);
        [u,s,v] = svdsim(S*(V'*U_new)*S_new);
        U = U*u;
        S = s;
        V = V_new*v;

    end
    %B_L2_inv = U * S * V';
end


