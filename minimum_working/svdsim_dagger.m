function [Q,D,R] = svdsim_dagger(A)
   %It is named svdsim only because it looks like svd and it works like svd in
% RDQ
    [Q_d,D_d,R_d] = svdsim(A');
     Q = R_d;
     D = D_d';
     R = Q_d;
end

