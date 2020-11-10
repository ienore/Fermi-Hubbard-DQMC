function B_L_inv = Get_B_L_inv(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    %Diag = eye([NumOfVertexs,NumOfVertexs])*(Miu-Uene/2);
    Sigma_L = Sigma(L,:);
    V_L = zeros([NumOfVertexs,NumOfVertexs]);
    for i = 1:1:NumOfVertexs
       V_L(i,i) = exp(-lambda*alpha*Sigma_L(i));
    end %Buide Up Matrix V_L
    B_L_inv = expm(-D_Tau*(K*T_hop))*V_L;
%     if L == 0 || L == TempSlice
%         B_L_inv = V_L;
%     end
end


