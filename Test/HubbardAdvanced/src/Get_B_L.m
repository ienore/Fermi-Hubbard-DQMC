function B_L = Get_B_L(alpha,L)
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
    V_L = diag(exp(lambda*alpha*Sigma(L,:)));
    B_L = V_L*expmK;
end




