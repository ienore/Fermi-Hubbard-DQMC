function G_L = Get_G_L_from_list(alpha,L)
    global TempSlice;
    global N_wrap;
    global u_up_list_1;
    global d_up_list_1;
    global v_up_list_1;
    global u_down_list_1;
    global d_down_list_1;
    global v_down_list_1;
    global u_up_list_2;
    global d_up_list_2;
    global v_up_list_2;
    global u_down_list_2;
    global d_down_list_2;
    global v_down_list_2;
    
    if mod(TempSlice,N_wrap) ~= 0
       fprintf("Get_G_L_from_list ERROR\n");
    end 

    refine_index = L / N_wrap;
    
    if alpha > 0
        U_R = u_up_list_1(:,:,refine_index);
        D_R = d_up_list_1(:,:,refine_index);
        V_R = v_up_list_1(:,:,refine_index);
        V_L = u_up_list_2(:,:,refine_index);
        D_L = d_up_list_2(:,:,refine_index);
        U_L = v_up_list_2(:,:,refine_index);
    else
        U_R = u_down_list_1(:,:,refine_index);
        D_R = d_down_list_1(:,:,refine_index);
        V_R = v_down_list_1(:,:,refine_index);
        V_L = u_down_list_2(:,:,refine_index);
        D_L = d_down_list_2(:,:,refine_index);
        U_L = v_down_list_2(:,:,refine_index);
    end

    V_R = V_R'; 
    U_L = U_L';

    [u,d,v] = svdsim(D_R*(V_R*V_L)*D_L);
    v = v';
    U = U_R*u;
    D = d;
    V = v*U_L;
    [u,d,v] = svdsim(U'/V+D);
    v = v';
    d_vec = diag(d);
    d_inv = diag(1.0./d_vec,0);

    G_L = (v*V)\d_inv/(U*u);

   
end


