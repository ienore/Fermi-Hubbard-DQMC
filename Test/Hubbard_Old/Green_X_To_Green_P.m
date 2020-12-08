function Green_P = Green_X_To_Green_P(Green_X,NumInEdge,NumOfVertexs)
% No factor 2, Please input only one spin direction Green Function
    for px_index = 1:1:NumInEdge
        for py_index = 1:1:NumInEdge
            delta_p = 2*pi/(NumInEdge);
            px = (px_index-NumInEdge) * delta_p + pi ;
            py = (py_index-NumInEdge) * delta_p + pi;
            temp_sum = 0.0;
            for site_m = 1:1:NumOfVertexs
                for site_n = 1:1:NumOfVertexs
                    [mx,my] = IndexToCoor(site_m,NumInEdge);
                    [nx,ny] = IndexToCoor(site_n,NumInEdge);
                    delta_x = mx - nx;
                    delta_y = my - ny;
                    temp_sum = temp_sum + exp(1i * (px*delta_x + py*delta_y))*Green_X(site_m,site_n);
                end
            end
            Green_P(px_index,py_index) = real(temp_sum) / (NumOfVertexs);%% Why ?
        end
    end
end

