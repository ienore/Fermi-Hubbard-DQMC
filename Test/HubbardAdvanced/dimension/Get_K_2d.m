function K = Get_K_2d(NumInEdge)
    NumOfVertexs = NumInEdge^2;
    K = zeros([NumOfVertexs,NumOfVertexs]);
    for i = 1:1:NumOfVertexs
        for j = 1:1:NumOfVertexs
            [xi,yi] = IndexToCoor_2d(i,NumInEdge);
            [xj,yj] = IndexToCoor_2d(j,NumInEdge);
            if ((min(abs(xi-xj),abs(abs(xi-xj)-NumInEdge)))+(min(abs(yi-yj),abs(abs(yi-yj)-NumInEdge)))) == 1
                K(i,j) = 1;
            end
        end
    end%Buide Up Matrix K
end

