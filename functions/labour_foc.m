function N=labour_foc(C,K,A,par)

num=(1-par.alpha)*A*C^(-par.gamma)*K^par.alpha;
deno=par.theta;
power=1/(par.varphi+par.alpha);

N=(num/deno)^power;

end 