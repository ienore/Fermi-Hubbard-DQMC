function [Q,D,R] = svdsim(A)
   %It is named svdsim only because it looks like svd and it works like svd in
%our DQMC algorithm, but in fact, this is qr based UDT decomposition algorithum
    N = length(A);
    [Q,R] = qr(A);
    D = zeros(N);
    for site_index = 1:1:N
        D(site_index,site_index) = R(site_index,site_index);
        R(site_index,:) = R(site_index,:)/D(site_index,site_index);
    end
    R = R';
end

