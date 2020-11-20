function B_L = Get_B_L(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    global expmK;
    %Diag = eye([NumOfVertexs,NumOfVertexs])*(Miu-Uene/2);
    V_L = diag(lambda*alpha*Sigma(L,:));
     B_L = V_L*expmK;
%     if L == 0 || L == TempSlice
%         B_L = V_L;
%     end
end


