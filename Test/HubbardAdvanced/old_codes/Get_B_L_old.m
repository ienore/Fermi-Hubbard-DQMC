function B_L = Get_B_L_old(alpha,L)
	global NumInEdge;
    global NumOfVertexs;
    global D_Tau;
    global TempSlice;
    global Uene;
    global Miu;
    global Beta;
    global lambda;
    global T_hop;
    global N_wrap;
    global id_mat;
    global Sigma;
    global expmK;
	global B_up_list;
	global B_down_list;
    %Diag = eye([NumOfVertexs,NumOfVertexs])*(Miu-Uene/2);
    Sigma_L = Sigma(L,:);
    V_L = zeros([NumOfVertexs,NumOfVertexs]);
    for i = 1:1:NumOfVertexs
       V_L(i,i) = exp(lambda*alpha*Sigma_L(i));
    end %Buide Up Matrix V_L
     B_L = V_L*expm(D_Tau*(K*T_hop));
%     if L == 0 || L == TempSlice
%         B_L = V_L;
%     end
end


