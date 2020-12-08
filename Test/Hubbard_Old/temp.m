t_start = 3;
propa_time_index = 8;

G_test = Get_G_L2(1,propa_time_index,t_start,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
G_equal = Get_G_L(1.0,t_start,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);%up state
B_propa = Get_B_L2(1,propa_time_index,t_start,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene); 
G_propa = B_propa * G_equal;

G_diff = G_test - G_propa;
disp(norm(G_diff))