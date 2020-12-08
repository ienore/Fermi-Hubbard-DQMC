function O_alpha = Get_O_alpha_p(alpha,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,TempSlice_0,p,K,T_hop)
    O_alpha = zeros([NumOfVertexs*p,NumOfVertexs*p]);
    temp_B = eye([NumOfVertexs,NumOfVertexs]);
    for p_index = 1:1:p
        temp_B = eye([NumOfVertexs,NumOfVertexs]);
        for L = (TempSlice_0*(p_index-1)+1):1:p_index*TempSlice_0
            Sigma_L = Sigma(L,:);
            V_L = zeros([NumOfVertexs,NumOfVertexs]);
            for i = 1:1:NumOfVertexs
               V_L(i,i) = Sigma_L(i);
            end %Buide Up Matrix V_L
            B_L = expm(-D_Tau*K*T_hop)*expm(alpha*lambda*V_L);
            temp_B = B_L*temp_B;
        end
        if (p_index == p)
            O_alpha(1:NumOfVertexs,(p_index-1)*NumOfVertexs+1:p_index*NumOfVertexs) = temp_B;
        else
            O_alpha(p_index*NumOfVertexs+1:(p_index+1)*NumOfVertexs,(p_index-1)*NumOfVertexs+1:p_index*NumOfVertexs) = -temp_B;
        end
    end
    O_alpha = O_alpha + eye([NumOfVertexs*p,NumOfVertexs*p]);
         
end


