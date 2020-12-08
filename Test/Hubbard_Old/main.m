% H = -t(c^\dagger c + H.C.) + U*n1*n2 - \miu (n1+n2) (From µÍÎÂËã·¨.pdf)

T_hop = -1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumInEdge = 4;
NumOfVertexs = NumInEdge^2;
%NumOfVertexs = 6;
K = Get_K(NumInEdge);
%K = Get_K_1d(NumOfVertexs);
%lambda = acosh(exp(Uene*D_Tau/2.0));


FINAL_RUN = 16;
NumOfWarm = 300;
NumOfEpoch = 2000;
U_range = 0:0.5:7;
Res_range = U_range*0.0;
Var_range = U_range*0.0;
Beta = 4.0;
for mi = 1:1:length(U_range)
    Uene = U_range(mi);
    Miu = Uene/2;
    D_Tau = 0.125;
    TempSlice = Beta/D_Tau;
    lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
    [return_mean,return_svar]=Sub_Par(mi,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);
    Res_range(mi) = return_mean;
    Var_range(mi) = return_svar;
    fprintf("0.025 - Uene = %f,  Mean = %f, Var_Ratio = %f\n",Uene,return_mean,abs(return_svar/return_mean));
end
errorbar(U_range,Res_range,Var_range);
xlabel("U/t");
ylabel("Structure Factor")
legend("D_Tau = 0.125")
hold on
    