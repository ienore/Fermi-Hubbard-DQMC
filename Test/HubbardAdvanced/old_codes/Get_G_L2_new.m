function G_L2 = Get_G_L2_new(alpha,L)
    global NumOfVertexs;
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
    
    refine_index = L / N_wrap;
    
    if alpha > 0
        Q_R = u_up_list_1(:,:,refine_index);
        D_R = d_up_list_1(:,:,refine_index);
        T_R = v_up_list_1(:,:,refine_index);
        T_L = u_up_list_2(:,:,refine_index);
        D_L = d_up_list_2(:,:,refine_index);
        Q_L = v_up_list_2(:,:,refine_index);
    else
        Q_R = u_down_list_1(:,:,refine_index);
        D_R = d_down_list_1(:,:,refine_index);
        T_R = v_down_list_1(:,:,refine_index);
        T_L = u_down_list_2(:,:,refine_index);
        D_L = d_down_list_2(:,:,refine_index);
        Q_L = v_down_list_2(:,:,refine_index);
    end
    
%     [Q_R,D_R,T_R] = Get_B_L2_svd_old(alpha,L,0);
%     [T_L,D_L,Q_L] = Get_B_L2_svd_old(alpha,TempSlice,L);
    T_R = T_R';
   Q_L = Q_L';
   D_R_min = zeros([NumOfVertexs]);
   D_R_max = zeros([NumOfVertexs]);
   D_L_min = zeros([NumOfVertexs]);
   D_L_max = zeros([NumOfVertexs]);
   D_R_max_inv = zeros([NumOfVertexs]);
   D_L_max_inv = zeros([NumOfVertexs]);

   for i = 1:1:NumOfVertexs
      if D_R(i,i) <= 1.0
          D_R_max(i,i) = 1.0;
          D_R_min(i,i) = D_R(i,i);
      else
          D_R_min(i,i) = 1.0;
          D_R_max(i,i) = D_R(i,i);
      end
      if D_L(i,i) <= 1.0
          D_L_max(i,i) = 1.0;
          D_L_min(i,i) = D_L(i,i);
      else
          D_L_min(i,i) = 1.0;
          D_L_max(i,i) = D_L(i,i);
      end
      D_R_max_inv(i,i) = 1.0/D_R_max(i,i);
      D_L_max_inv(i,i) = 1.0/D_L_max(i,i);

   end
   temp = D_R_max_inv /(Q_L*Q_R)*D_L_max_inv + D_R_min*(T_R*T_L)*D_L_min;
   G_L2 = Q_L\D_L_max_inv/temp*D_R_min*T_R;

    

end


