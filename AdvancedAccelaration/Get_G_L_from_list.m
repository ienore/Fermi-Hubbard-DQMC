function G_L = Get_G_L_from_list(alpha,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene)
    global N_wrap;
    global B_up_list;
    global B_down_list;
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
    global id_mat;
    refine_index_low = fix(L/N_wrap);
    refine_index_high = refine_index_low + 1;
    Binned_B_1 = eye(NumOfVertexs);
    Binned_B_2 = eye(NumOfVertexs);
    
    %% Best Situation:����
    if mod(L,N_wrap) == 0
        U_R = u_up_list_1(:,:,refine_index_low);
        D_R = d_up_list_1(:,:,refine_index_low);
        V_R = v_up_list_1(:,:,refine_index_low);
        U_L = u_up_list_2(:,:,refine_index_low);
        D_L = d_up_list_2(:,:,refine_index_low);
        V_L = v_up_list_2(:,:,refine_index_low);
    else
       %% Common Situation
        if refine_index_high * N_wrap < TempSlice && refine_index_low > 0
            if alpha > 0
                u_R = u_up_list_1(:,:,refine_index_low);
                d_R = d_up_list_1(:,:,refine_index_low);
                v_R = v_up_list_1(:,:,refine_index_low);
                u_L = u_up_list_2(:,:,refine_index_high);
                d_L = d_up_list_2(:,:,refine_index_high);
                v_L = v_up_list_2(:,:,refine_index_high);

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_up_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:N_wrap * refine_index_high
                    Binned_B_2 = B_up_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            else
                u_R = u_down_list_1(:,:,refine_index_low);
                d_R = d_down_list_1(:,:,refine_index_low);
                v_R = v_down_list_1(:,:,refine_index_low);
                u_L = u_down_list_2(:,:,refine_index_high);
                d_L = d_down_list_2(:,:,refine_index_high);
                v_L = v_down_list_2(:,:,refine_index_high);

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_down_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:N_wrap * refine_index_high
                    Binned_B_2 = B_down_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            end
        end
       %% if being in the highest bin
        if refine_index_high * N_wrap >= TempSlice
            if alpha > 0
                u_R = u_up_list_1(:,:,refine_index_low);
                d_R = d_up_list_1(:,:,refine_index_low);
                v_R = v_up_list_1(:,:,refine_index_low);
                u_L = id_mat;
                d_L = id_mat;
                v_L = id_mat;

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_up_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:TempSlice
                    Binned_B_2 = B_up_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            else
                u_R = u_down_list_1(:,:,refine_index_low);
                d_R = d_down_list_1(:,:,refine_index_low);
                v_R = v_down_list_1(:,:,refine_index_low);
                u_L = id_mat;
                d_L = id_mat;
                v_L = id_mat;

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_down_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:TempSlice
                    Binned_B_2 = B_down_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            end
        end
        %% if being in the lowest bin
        if refine_index_low == 0
            if alpha > 0
                u_R = id_mat;
                d_R = id_mat;
                v_R = id_mat;
                u_L = u_up_list_2(:,:,refine_index_high);
                d_L = d_up_list_2(:,:,refine_index_high);
                v_L = v_up_list_2(:,:,refine_index_high);

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_up_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:N_wrap * refine_index_high
                    Binned_B_2 = B_up_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            else
                u_R = id_mat;
                d_R = id_mat;
                v_R = id_mat;
                u_L = u_down_list_2(:,:,refine_index_high);
                d_L = d_down_list_2(:,:,refine_index_high);
                v_L = v_down_list_2(:,:,refine_index_high);

                for index = (refine_index_low * N_wrap+1):1:L
                    Binned_B_1 = B_down_list(:,:,index) * Binned_B_1;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_1);
                [u,s,v] = svdsim(S_new*(V_new'*u_R)*d_R);
                U_R = U_new*u;
                D_R = s;
                V_R = v_R*v;

                for index = (L+1):1:N_wrap * refine_index_high
                    Binned_B_2 = B_down_list(:,:,index) * Binned_B_2;
                end
                [U_new,S_new,V_new] = svdsim(Binned_B_2);
                [u,s,v] = svdsim(d_L*(v_L'*U_new)*S_new);
                U_L = u_L * u;
                D_L = s;
                V_L = V_new * v;
            end
        end
    end
    
    
    
  
   %[U_R,D_R,V_R] = Get_B_L2_svd(alpha,L,0,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
   %[V_L,D_L,U_L] = Get_B_L2_svd_dagger(alpha,TempSlice,L,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
  % G_L = inv(eye(NumOfVertexs) + U_R*D_R*V_R'*V_L*D_L*U_L');
   % I + U_R( D_R V_R V_L D_L) U_L = I + U_R u d v U_L
   V_R = V_R'; 
   U_L = U_L';
   
   [u,d,v] = svdsim(D_R*(V_R*V_L)*D_L);
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
   G_L = U_new*D_new*V_new';
end
