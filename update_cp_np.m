function [cprime, nprime] = update_cp_np(policy_k, k_grid_target, kprime_grid, grillA, par)
    %Guess for c_prime and n_prime
    nprime = 0.3 * ones(par.N, par.nkap);
    cprime = zeros(par.N, par.nkap);    
    for i = 1:par.N
        A = grillA(i);
        for j = 1:par.nkap
            k_p = kprime_grid(j);
            k_pp = interp1(k_grid_target, policy_k(i, :), k_p,'linear', 'extrap');
            n_old = 1/3;
            for it = 1:200
                y = A * k_p^par.alpha * n_old^(1-par.alpha);
                c_new = y + (1-par.delta)*k_p - k_pp;
                c_new = max(c_new, 1e-10);
                n_new = ((1-par.alpha)*A*k_p^par.alpha*(c_new)^(-par.gamma))^(1/(par.gamma+par.alpha));
                if abs(n_new - n_old) < 1e-6
                    break;
                end
                n_old = n_new;
            end
            nprime(i,j) = n_new;
            cprime(i,j) = c_new;
        end
    end
end
