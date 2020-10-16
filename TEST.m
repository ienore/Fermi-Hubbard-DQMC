N = 8;
[Q,R] = qr(A);
D = zeros([N]);
for site_index = 1:1:N
    D(site_index,site_index) = R(site_index,site_index);
    R(site_index,:) = R(site_index,:)/D(site_index,site_index);
end
R = R';
A_qr = Q*D*R';
diff = A_qr - A;
