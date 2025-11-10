function k_prime=lom_function(I,K, par)

part=(par.phi/2)*(I/K-par.delta)^2*K;

k_prime=I-part+(1-par.delta)*K;

end 