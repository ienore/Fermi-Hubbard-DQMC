function K = Get_K_3d(NumInEdge)
    NumOfVertexs = NumInEdge^3;
    K = zeros([NumOfVertexs,NumOfVertexs]);
    for i = 1:1:NumOfVertexs
        for j = 1:1:NumOfVertexs
            [xi,yi,zi] = IndexToCoor_3d(i,NumInEdge);
            [xj,yj,zj] = IndexToCoor_3d(j,NumInEdge);
            delta_x = min(abs(xi-xj),NumInEdge - abs(xi-xj));
            delta_y = min(abs(yi-yj),NumInEdge - abs(yi-yj));
            delta_z = min(abs(zi-zj),NumInEdge - abs(zi-zj));
            if (delta_x + delta_y + delta_z) == 1
                K(i,j) = 1;
            end
        end
    end%Buide Up Matrix K
end

