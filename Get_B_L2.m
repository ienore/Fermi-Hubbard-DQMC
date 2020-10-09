function B_L2 = Get_B_L2(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2 = eye(NumOfVertexs);
    for L = L1+1:1:L2
        B_L2 = Get_B_L(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*B_L2;
    end
end


