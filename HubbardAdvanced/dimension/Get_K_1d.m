function K = Get_K_1d(NumOfVertexs)
    K = zeros([NumOfVertexs,NumOfVertexs]);
    for i = 1:1:NumOfVertexs
       if i==1
           K(i,i+1) = 1;
           K(i,NumOfVertexs) = 1;
       elseif i == NumOfVertexs
           K(i,i-1) = 1;
           K(i,1) = 1;
       else 
           K(i,i-1) = 1;
           K(i,i+1) = 1;
       end
    end
end

