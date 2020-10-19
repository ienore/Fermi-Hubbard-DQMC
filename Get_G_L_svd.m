function [U_new,D_new,V_new] = Get_G_L_svd(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
%       G_L = inv(eye(NumOfVertexs) + Get_B_L2(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)* ...
%                                     Get_B_L2(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene));
   
   [U_R,D_R,V_R] = Get_B_L2_svd(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   [V_L,D_L,U_L] = Get_B_L2_svd_dagger(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
  % G_L = inv(eye(NumOfVertexs) + U_R*D_R*V_R'*V_L*D_L*U_L');
   % I + U_R( D_R V_R V_L D_L) U_L = I + U_R u d v U_L
   V_R = V_R'; 
   U_L = U_L';
   
   [u,d,v] = svdsim(D_R*V_R*V_L*D_L);
   v = v';
   U = U_R*u;
   D = d;
   V = v*U_L;
   [u,d,v] = svdsim(U'/V+D);
   v = v';
    d_vec = diag(d);
    d_inv = diag(1.0./d_vec,0);
   
    U_new = inv(v*V);
    D_new = d_inv;
    V_new = inv(U*u)';
   %G_L = U_new*D_new*V_new';
end