function B_L2_inv = Get_B_L2_inv(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2_inv = eye(NumOfVertexs);
    for L = L1+1:1:L2
        B_L2_inv = B_L2_inv*Get_B_L_inv(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
    end
end


