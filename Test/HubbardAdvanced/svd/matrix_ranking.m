function [A_P,P] = matrix_ranking(A)
    N = length(A);
    norm_value = zeros([1,N]);
    for site_index = 1:1:N
       norm_value(site_index) = norm(A(:,site_index)); 
    end
    [sorted,idx] = sort(-norm_value);
    P = zeros(N);
    for site_index = 1:1:N
        P(idx(site_index),site_index) = 1;
    end
    A_P = A*P;
end

