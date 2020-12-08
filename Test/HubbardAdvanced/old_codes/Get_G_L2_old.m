function G_L2 = Get_G_L2_old(alpha,L)
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
    [Q_R,D_R,T_R] = Get_B_L2_svd_old(alpha,L,0);
    [T_L,D_L,Q_L] = Get_B_L2_svd_old(alpha,TempSlice,L);
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


