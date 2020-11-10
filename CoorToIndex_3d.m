function index = CoorToIndex_3d(x,y,z,NumInEdge)
    index = NumInEdge^2*(z-1) + NumInEdge*(y-1)+x;
end

