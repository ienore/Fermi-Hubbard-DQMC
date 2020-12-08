function Show = ReArrangeToShow(NumInEdge,Input)
    for i = 1:1:NumInEdge^2
        xi_n = fix(double(i-1)/NumInEdge)+1;
        yi_n = mod(double(i-1),NumInEdge)+1;
        Show(xi_n,yi_n)=Input(i);
    end
end

