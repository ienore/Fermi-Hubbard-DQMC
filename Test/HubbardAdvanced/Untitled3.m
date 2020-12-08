Uene = 4;
Energy = [0.1,-0.4,-0.85,-1.05];
E_k = [-0.3336,-0.54014,-0.72189,-0.82583];
DD = Uene*[0.1935,0.165848,0.155357,0.158064];
N = 1e6;
step = 0.01;
a = 0;
b = 0;
c = 0;
diffe_old = 10^3;
for index = 1:1:N
    if mod(index,N/100) == 0
        fprintf("ratio = %f\t diff = %e\t a=%f,b=%f,c=%f\n",index/N,diffe_old,a,b,c);
    end
    a_g = a + step * (rand-0.5);
    b_g = b + step * (rand-0.5);
    c_g = c + step * (rand-0.5);
    Energy_guess = a_g * E_k + b_g * DD + c_g*Uene;
    diffe = norm(Energy_guess - Energy);
    diffe_old = norm(a* E_k + b * DD + c*Uene - Energy);
    if diffe < diffe_old
        a = a_g;
        b = b_g;
        c = c_g;
    end
end
    
2* E_k + 1 * DD + 0*Uene