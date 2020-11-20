function [x,y] = IndexToCoor(index,NumInEdge)
    x = fix(double(index-1)/NumInEdge)+1;
    y = mod(double(index-1),NumInEdge)+1;
end

