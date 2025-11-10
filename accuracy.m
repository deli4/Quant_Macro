function ee = accuracy(par, grids, c_pol, n_pol, k_pol, transition, marg_u_c)
% grids.k : capital grid (strictly increasing)
% grids.a : discretized productivity states (5x1). Use exp() below if logs.

state = size(transition,2);   % assumes square P
N     = numel(grids.k);
ee    = nan(N, state);

for a = 1:state
  for k = 1:N
    Ccur = c_pol(a,k);
    if Ccur <= 0, ee(a,k) = NaN; 
        continue; 
    end

    muc = marg_u_c(Ccur);
    Kp  = k_pol(a,k);
    if Kp <= grids.k(1) || Kp >= grids.k(end), 
        ee(a,k) = NaN; 
        continue; 
    end

    rhs = 0.0;
    for ap = 1:state
      Cn = linear_interpolation(grids.k, c_pol(ap,:), Kp);
      Nn = linear_interpolation(grids.k, n_pol(ap,:), Kp);
      if Cn <= 0 || Nn <= 0 || Nn > 1, 
          continue; 
      end

      % If grids.a are logs, use: Aap = exp(grids.a(ap));
      Aap = grids.a(ap);

      Rn  = 1 - par.delta + par.alpha * Aap * (Kp^(par.alpha-1)) * (Nn^(1-par.alpha));
      rhs = rhs + transition(ap,a) * par.beta * (marg_u_c(Cn)/muc) * Rn;
    end

    ee(a,k) = abs(1 - rhs);
  end
end

% Plot in log10 scale by productivity state
figure; hold on;
for a = 1:state
  plot(grids.k, log10(ee(a,:) + eps), 'DisplayName', sprintf('A_%d', a));
end
hold off;
xlabel('Capital k');
ylabel('log_{10}(|Euler residual|)');
title('Euler equation errors by productivity state, \phi=', par.phi);
legend('show','Location','best'); grid on;
end
