function O_alpha = Get_O_alpha(alpha,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    O_alpha = zeros([NumOfVertexs*TempSlice,NumOfVertexs*TempSlice]);
    Diag = eye([NumOfVertexs,NumOfVertexs])*(Miu-Uene/2);
    for L = 1:1:TempSlice
        Sigma_L = Sigma(L,:);
        V_L = zeros([NumOfVertexs,NumOfVertexs]);
        for i = 1:1:NumOfVertexs
           V_L(i,i) = lambda*alpha*Sigma_L(i);
        end %Buide Up Matrix V_L
        B_L = expm(-D_Tau*(-K*T_hop-Diag))*expm(V_L);
        if(L == TempSlice)
            O_alpha(1:NumOfVertexs,(TempSlice-1)*NumOfVertexs+1:NumOfVertexs*TempSlice) = B_L;
        else
            O_alpha(L*NumOfVertexs+1:(L+1)*NumOfVertexs,(L-1)*NumOfVertexs+1:L*NumOfVertexs) = - B_L;
        end
    end
    O_alpha = O_alpha + eye([NumOfVertexs*TempSlice,NumOfVertexs*TempSlice]);
end


