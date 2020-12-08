function B_L_inv = Get_B_L_inv(alpha,L)
    global lambda;
    global Sigma;    
    global expmK_inv;
    V_L = diag(exp(-lambda*alpha*Sigma(L,:)));
    B_L_inv = expmK_inv*V_L;
end


