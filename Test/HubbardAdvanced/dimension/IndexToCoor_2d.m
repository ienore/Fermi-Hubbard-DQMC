function [x,y] = IndexToCoor_2d(index,NumInEdge)
    x = fix(double(index-1)/NumInEdge)+1;
    y = mod(double(index-1),NumInEdge)+1;
end

