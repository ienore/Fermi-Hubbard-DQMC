function B_L_inv = Get_B_L_inv(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    %Diag = eye([NumOfVertexs,NumOfVertexs])*(Miu-Uene/2);
    global expmK_inv;
    V_L = diag(-lambda*alpha*Sigma(L,:));
    B_L_inv = expmK_inv*V_L;
%     if L == 0 || L == TempSlice
%         B_L_inv = V_L;
%     end
end


