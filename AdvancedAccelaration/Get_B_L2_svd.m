function [U,S,V] = Get_B_L2_svd(alpha,L2,L1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    B_L2 = eye(NumOfVertexs);
    [U,S,V] = svdsim(B_L2);
    Bin_Size = 7;
    global B_up_list;
    global B_down_list;
    
    if fix((L2-L1)/Bin_Size)==(L2-L1)/Bin_Size
        Bin_run = (L2-L1)/Bin_Size;
    else
        Bin_run = fix((L2-L1)/Bin_Size)+1;
    end
    
    for Bin_index = 1:1:Bin_run
        Binned_B_L = eye(NumOfVertexs);
        L_start = L1 + (Bin_index-1) * Bin_Size ;
        if Bin_index < Bin_run
            for index = 1:1:Bin_Size
                if alpha > 0
                    B = B_up_list(:,:,L_start+index);
                else
                    B = B_down_list(:,:,L_start+index);
                end
                Binned_B_L = B * Binned_B_L;
                %Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
            end
        else
            for index = 1:1:(L2 - L_start)
                if alpha > 0
                    B = B_up_list(:,:,L_start+index);
                else
                    B = B_down_list(:,:,L_start+index);
                end
                Binned_B_L = B * Binned_B_L;
                %Binned_B_L = Get_B_L(alpha,L_start + index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)*Binned_B_L;
            end
        end
        
        [U_new,S_new,V_new] = svdsim(Binned_B_L);
        [u,s,v] = svdsim(S_new*(V_new'*U)*S);
        U = U_new*u;
        S = s;
        V = V*v;
       % [U,S,V_new] = svdsim(Binned_B_L*U*S);
       % V = V * V_new;
    end
    %B_L2 = U * S * V';
end

