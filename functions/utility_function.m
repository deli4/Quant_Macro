function u=utility_function(C,N,par)

u_c=(C^(1-par.gamma)-1)/(1-par.gamma);
u_n=(par.theta/(1+par.varphi))*N^(1+par.varphi);

u=u_c-u_n;

end 