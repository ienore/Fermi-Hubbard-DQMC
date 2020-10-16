alpha = 1;
L1 = 0;




%fprintf("FINAL_RUN:%d\n",FINAL_RUN);

count = 1;
T_hop = 1.0;
NumInEdge = 7;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

Uene = 6;
Miu = Uene/2;
Beta = 6;
D_Tau = 0.2;
TempSlice = Beta/D_Tau;
L2 = TempSlice;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
sample = zeros([NumOfVertexs,TempSlice]);

B_L2 = eye(NumOfVertexs);
[U,S,V] = svdsim(B_L2);
Bin_Size = 7;
NumBin = fix((L2-L1)/Bin_Size)+1;
if mod(L2-L1,Bin_Size)==0
    NumBin = (L2-L1)/Bin_Size;
end
L2 = TempSlice;
for Bin_index = 1:1:NumBin
    Binned_B_L = eye(NumOfVertexs);
    L_start = L1 + (Bin_index-1) * Bin_Size ;
    if Bin_index < NumBin    
        for index = 1:1:Bin_Size
            Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
        end
    else
        for index = 1:1:(L2 - L_start)
            Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
        end
    end
    [U_t,S_t,V_t] = svdsim(Binned_B_L);
    [U_n,S_n,V_n] = svdsim(S_t*(V_t'*U)*S);
    U = U_t*U_n;
    V = V * V_n;
    S = S_n;
    [U,S,V_new] = svdsim((Binned_B_L*U)*S);
    V = V * V_new;
    for site_index = 1:1:NumOfVertexs
        sample(site_index,Bin_index) = S(site_index,site_index);
    end
     if Bin_index ==4
         break;
     end
end
B_direct = (Binned_B_L*U)*S;
[U_a,S_a,V_a] =svdsim((Binned_B_L*U)*S);
B_svdsim_ed = U_a*S_a*V_a';
B_diff = B_direct - B_svdsim_ed;
B_L2 = U * S * V';

% for index = L1+1:1:L2
%     B_L2 = Get_B_L(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*B_L2;
% end

% [U,S,V] = svdsim(B_L2);