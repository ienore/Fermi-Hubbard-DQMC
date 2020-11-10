NumInEdge = 5;
for x = 1:1:NumInEdge
    for y = 1:1:NumInEdge
        for z = 1:1:NumInEdge
            index = CoorToIndex_3d(x,y,z,NumInEdge);
            [xp,yp,zp] = IndexToCoor_3d(index,NumInEdge);
            fprintf("%d,%d,%d---%d,%d,%d\n",x,y,z,xp,yp,zp);
        end
    end
end