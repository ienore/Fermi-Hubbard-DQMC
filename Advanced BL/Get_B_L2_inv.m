function B_L2_inv = Get_B_L2_inv(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2_inv = eye(NumOfVertexs);
    [U,S,V] = svd(B_L2_inv);
    Bin_Size = 5;
    for Bin_index = 1:1:fix((L2-L1)/Bin_Size)+1
        Binned_B_L = eye(NumOfVertexs);
        L_start = L1 + (Bin_index-1) * Bin_Size ;
        if Bin_index < fix((L2-L1)/Bin_Size)+1       
            for index = 1:1:Bin_Size
                Binned_B_L = Binned_B_L*Get_B_L_inv(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            end
        else
            for index = 1:1:(L2 - L_start)
                Binned_B_L = Binned_B_L*Get_B_L_inv(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            end
        end
        [U_new,S,V] = svd(S * (V'* Binned_B_L));
        U = U*U_new;
    end
    B_L2_inv = U * S * V';
end


