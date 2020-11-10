Tau_range = 0:D_Tau:Beta;
px = pi;
py = pi;
Miu = 0;
e_k = - 2*T_hop*(cos(px)+cos(py));
pre_factor = 1.0-1.0/(exp(Beta*(e_k-Miu))+1);
plot_theo = pre_factor * exp(-e_k*Tau_range);
plot_1 = plot_theo;

px = 0;
py = 0;
Miu = 0;
e_k = - 2*T_hop*(cos(px)+cos(py));
pre_factor = 1.0-1.0/(exp(Beta*(e_k-Miu))+1);
plot_theo = pre_factor * exp(-e_k*Tau_range);
plot_2 = plot_theo;


plot(Tau_range,log(plot_1 + plot_2),'blue');

    