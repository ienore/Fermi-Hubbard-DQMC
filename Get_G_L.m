function G_L = Get_G_L(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
%       G_L = inv(eye(NumOfVertexs) + Get_B_L2(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)* ...
%                                     Get_B_L2(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene));
   [U_R,D_R,V_R] = Get_B_L2_svd(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   [V_L,D_L,U_L] = Get_B_L2_svd(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   
   V_R = V_R';
   U_L = U_L';
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
   
   temp = D_R_max_inv/(U_L*U_R)*D_L_max_inv+D_R_min*V_R*V_L*D_L_min;
   G_L = U_L\D_L_max_inv/temp*D_R_max_inv/U_R;
   
   
   %G_L = inv(U_L)*D_L_max_inv*inv(D_R_max_inv*inv(U_L*U_R)* ...
   %    D_L_max_inv + D_R_min*V_R*V_L + ...
   %    D_L_min)*D_R_max_inv*inv(U_R);%%This method is more unstable

   
end


