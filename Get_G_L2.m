function G_L2 = Get_G_L2(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
%      G_L = inv(eye(NumOfVertexs) + Get_B_L2(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)* ...
%                                    Get_B_L2(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene));
   B = Get_B_L2(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   G_L = Get_G_L(alpha,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   G_L2 = B * G_L;%In PPt , this is the G_t

end


