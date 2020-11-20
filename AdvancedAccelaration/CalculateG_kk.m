NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);
T_hop = -1;
Uene = 8;
Miu = Uene/2;
Beta = 6;
D_Tau = 0.25;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
NumOfWarm = 50;
NumOfEpoch = 100;
zjy_index = 1;

Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 10;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_green_up_mat = zeros([NumOfVertexs,NumOfVertexs,(TempSlice-1)*NumOfEpoch*2]);
mea_result_green_down_mat = zeros([NumOfVertexs,NumOfVertexs,(TempSlice-1)*NumOfEpoch*2]);
mea_result_green_timely = zeros([NumOfVertexs,NumOfVertexs,(TempSlice-1)*NumOfEpoch*2]);
count = 1.0;
for warm_index = 1:1:NumOfWarm
    if mod(zjy_index,8) == 1
        fprintf("D_Tau = %f, Warming_Ratio = %f\n",D_Tau,warm_index/NumOfWarm);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index_mother,N_wrap) == 2
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
            end
        else
            time_index = TempSlice - time_index_mother + 1;
            if mod(time_index_mother,N_wrap) == 1
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
            end
        end
        for site_index = 1:1:NumOfVertexs
            delta_up = zeros(NumOfVertexs);
            delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
            delta_down = zeros(NumOfVertexs);
            delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
            R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
            R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
            PosToChange = R_up*R_down;
            if rand < PosToChange
                Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
            end
        end
    end
    end
end
for epoch_index = 1:1:NumOfEpoch
    if mod(zjy_index,8) == 1
        fprintf("D_Tau = %f,MC_Ratio = %f\n",D_Tau,epoch_index/NumOfEpoch);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index_mother,N_wrap) == 2
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
            end
        else
            time_index = TempSlice - time_index_mother + 1;
            if mod(time_index_mother,N_wrap) == 1
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        green_up = green_L_up;
        green_down = green_L_down;
        green_up_c = id_mat - transpose(green_L_up);%有没有可能本来的green就是conjugate过的
        green_down_c = id_mat - transpose(green_L_down);
        %mea_result_green_up_mat(:,:,count) = green_up_c;
        %mea_result_green_down_mat(:,:,count) = green_down_c;
        mea_result_green_timely(:,:,count) = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
        mea_result_green_up_mat(:,:,count) = green_up;
        mea_result_green_down_mat(:,:,count) = green_down;
        sum_measure = 0.0;
        
        mea_result(count) = sum_measure;
        mea_result_auxi(count) = green_down_c(1,1)*green_up_c(1,1);
        count = count + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for site_index = 1:1:NumOfVertexs
            delta_up = zeros(NumOfVertexs);
            delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
            delta_down = zeros(NumOfVertexs);
            delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
            R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
            R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
            PosToChange = R_up*R_down;
            if rand < PosToChange
                Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
            end
        end
    end
    end
end



mea_result = mea_result(1:(count-1));
mea_result_green_up_mat = mea_result_green_up_mat(:,:,1:(count-1));
mea_result_green_down_mat = mea_result_green_down_mat(:,:,1:(count-1));
NumInCut = fix(length(mea_result)/N_cut);
result_plot = zeros([1,N_cut]);
for cut_index = 1:1:N_cut
    result_plot(cut_index) = mean(mea_result(((cut_index-1)*NumInCut+1):cut_index*NumInCut));
end
return_mean = mean(mea_result);
return_svar = sqrt(var(result_plot));
%fprintf("S_mean = %f, S_std_var = %f\n",return_mean,return_svar)
plot_result_green_up_c = zeros(NumOfVertexs);
plot_result_green_down_c = zeros(NumOfVertexs);
for x_index = 1:1:NumOfVertexs
    for y_index = 1:1:NumOfVertexs
        plot_result_green_up_c(x_index,y_index) = mean(mea_result_green_up_mat(x_index,y_index,:));
        plot_result_green_down_c(x_index,y_index) = mean(mea_result_green_down_mat(x_index,y_index,:));
    end
end


green_up_c = plot_result_green_up_c;
green_down_c = plot_result_green_down_c;
green_c = green_up_c + green_down_c;
green_p_c = zeros(NumInEdge+1);
for px_index = 1:1:NumInEdge+1
    for py_index = 1:1:NumInEdge+1
        delta_p = 2*pi/(NumInEdge);
        px = (px_index-NumInEdge-1) * delta_p - pi ;
        py = (py_index-NumInEdge-1) * delta_p - pi;
        temp_sum = 0.0;
        for site_m = 1:1:NumOfVertexs
            for site_n = 1:1:NumOfVertexs
                [mx,my] = IndexToCoor(site_m,NumInEdge);
                [nx,ny] = IndexToCoor(site_n,NumInEdge);
                delta_x = mx - nx;
                delta_y = my - ny;
                temp_sum = temp_sum + exp(1i * (px*delta_x + py*delta_y))*green_c(site_m,site_n);
            end
        end
        green_p_c(px_index,py_index) = real(temp_sum) / (2*NumOfVertexs);%% Why ?
    end
end

[X,Y] = meshgrid(-pi:delta_p:pi);
[Xq,Yq] = meshgrid(-pi:0.1:pi);         
green_p_c_plot = interp2(X,Y,green_p_c,Xq,Yq,'cubic');
%surf(Xq,Yq,green_p_c_plot);
surf(X,Y,green_p_c);
xlabel("K_x");
ylabel("K_y");
zlabel("n(k)");
zlim([0,1]);
title(['Uene = ',num2str(Uene),'  ','Beta = ',num2str(Beta)])
colorbar;
caxis([0,1]);

% s = pcolor(green_p_c);
% xlabel("K_x");
% xlim([-pi,pi]);
% ylim([-pi,pi]);
% ylabel("K_y");
% title(['n(k) ','Uene = ',num2str(Uene),'  ','Beta = ',num2str(Beta)])
% s.FaceColor = 'interp';
% axis equal;
% colorbar;
% caxis([0,2])
% 

disp(mean(mea_result_auxi))

%[grad_x,grad_y] = gradient(green_p_c_plot);
%grad_sum = sqrt(grad_x.^2 .* grad_y.^2);
%surf(grad_sum);
