function G_L2 = Get_G_L2(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
%      G_L = inv(eye(NumOfVertexs) + Get_B_L2(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)* ...
%                                    Get_B_L2(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene));

%    B = Get_B_L2(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
%    G_L = Get_G_L(alpha,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
%    G_L2 = G*B;%In PPt , this is the G_t
%       [U_R,D_R,V_R] = Get_B_L2_svd(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
%       [U_L,D_L,V_L]= Get_G_L_svd(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);

%      [u,d,v] = svdsim(D_L*(V_L'*U_R)*D_R);
%      U_new = U_L*u;
%      D_new = d;
%      V_new = V_R*v;
%      G_L2 = U_new*D_new*V_new';
% 
%     [u,d,v] = svdsim((V_L'*U_R));
%      U_new = U_L*D_L*u;
%      D_new = d;
%      V_new = V_R*v*D_R';
%      G_L2 = U_new*D_new*V_new';
%      %G_L2 = U_L*(D_L*(V_L'*U_R)*D_R)*V_R';
%      %G_L2 = U_L*D_L*V_L'*U_R*D_R*V_R';
%      %G_L2 = G*B;
% 
% % 
    [Q_R,D_R,T_R] = Get_B_L2_svd(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
    [T_L,D_L,Q_L] = Get_B_L2_svd_dagger(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
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


