function [U,S,V] = Get_B_L2_svd(alpha,L2,L1)
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
	B_L2 = eye(NumOfVertexs);
    [U,S,V] = svdsim(B_L2);
    Bin_Size = N_wrap;
    
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
                Binned_B_L = B*Binned_B_L;
            end
        else
            for index = 1:1:(L2 - L_start)
				if alpha > 0 
					B = B_up_list(:,:,L_start+index);
				else
					B = B_down_list(:,:,L_start+index);
				end
                Binned_B_L = B*Binned_B_L;
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

