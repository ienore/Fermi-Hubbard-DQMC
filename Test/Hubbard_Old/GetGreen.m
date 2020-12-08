function Green = GetGreen(alpha,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop)
    O_alpha = eye([NumOfVertexs,NumOfVertexs]);
    
    for L = TempSlice:-1:1
        Sigma_L = Sigma(L,:);
        V_L = zeros([NumOfVertexs,NumOfVertexs]);
        for i = 1:1:NumOfVertexs
           V_L(i,i) = Sigma_L(i);
        end %Buide Up Matrix V_L
        B_L = expm(D_Tau*K*T_hop)*expm(alpha*lambda*V_L);
        O_alpha = O_alpha * B_L;
    end
    O_alpha = O_alpha + eye([NumOfVertexs,NumOfVertexs]);
    Green = inv(O_alpha);
end


