function OB = Get_OB(alpha,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop)
    global Uene;
    global Miu;
    temp = eye([NumOfVertexs,NumOfVertexs]);
    for L = TempSlice:-1:1
        Sigma_L = Sigma(L,:);
        V_L = zeros([NumOfVertexs,NumOfVertexs]);
        for i = 1:1:NumOfVertexs
           V_L(i,i) = lambda*alpha*Sigma_L(i)-D_Tau*(Miu-Uene/2.0);
        end %Buide Up Matrix V_L
        temp = temp*expm(-D_Tau*K*T_hop)*expm(V_L);
    end
    temp = temp + eye([NumOfVertexs,NumOfVertexs]);
    OB = inv(temp);
end


