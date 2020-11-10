function [U,S,V] = JOS_SVD(G)
%JOS_SVD  Jacobi One Sided SVD program
%
% A simple program that demonstrates how to use the jacobi method to obtain
% SVD
% A may be rectangular and complex.
%
% usage: [U,S,V]= JOS_SVD(A)
%     or      S = JOS_SVD(A)
%
% with A = U*S*V' , S>=0 , U'*U = Iu  , and V'*V = Iv
%
%see also: SVD, EIG, QR, BIDIAG, HESS
%

% Lateef Adewale Kareem
% January 7, 2018
error = 1;
n = size(G, 2);
V = eye(n); S = eye(n);
while error > 1e-6
    error = 0;
    for i = 1:n-1
        for j = i+1:n
            x = G(:,i);
            y = G(:,j);
            a = x'*x;
            b = y'*y;
            c = x'*y;
            
            error = max(error,abs(c)/sqrt(a*b));
            
            %Jacobi rotation
            e = (b - a)/(2*c);
            t = sign(e)/(abs(e) + sqrt(1 + e^2));
            cs = 1/sqrt(1+t^2);
            sn = cs*t;
            
            % update G
            x = G(:,i);
            G(:,i) = cs*x - sn*G(:,j);
            G(:,j) = sn*x + cs*G(:,j);
            
            % update V
            x = V(:,i);
            V(:,i) = cs*x - sn*V(:,j);
            V(:,j) = sn*x + cs*V(:,j);
            
        end
    end
end

U = G;
for i = 1:n
    S(i,i) = norm(G(:,i));
    U(:,i) = G(:,i)/S(i,i);
end

if nargout<=1
   U = S;
end