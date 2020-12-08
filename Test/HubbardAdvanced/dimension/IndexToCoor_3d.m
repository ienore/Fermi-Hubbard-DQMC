function [x,y,z] = IndexToCoor_3d(index,NumInEdge)
    x = mod(index-1,NumInEdge)+1;
    y = mod((index-x)/NumInEdge,NumInEdge)+1;
    z = (index - x - NumInEdge*(y-1))/NumInEdge^2+1;
end

