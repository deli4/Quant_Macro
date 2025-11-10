function [sim,m,sd,logcorr]=simulation_test(par,grid,P,ss,c_pol,n_pol,i_pol,k_pol,T,seed,shock_size,T_graph,outdir)

% 0) Plotting horizon guard
if nargin < 13 || isempty(T_graph)
    T_graph = min(40,T);
else
    T_graph = min(T_graph, T);
end

rng(seed);

% 1) Prealloc
sim.A=zeros(T,1); sim.C=zeros(T,1); sim.N=zeros(T,1);
sim.I=zeros(T,1); sim.Y=zeros(T,1); sim.K=zeros(T,1);
A_index=zeros(T,1);

% 2) Policy interpolants (linear in A, no K-extrap)
F_c = griddedInterpolant({grid.a, grid.k}, c_pol, 'linear','nearest');
F_n = griddedInterpolant({grid.a, grid.k}, n_pol, 'linear','nearest');
F_i = griddedInterpolant({grid.a, grid.k}, i_pol, 'linear','nearest');
F_k = griddedInterpolant({grid.a, grid.k}, k_pol, 'linear','nearest');

% 3) Initial conditions (one-time level shock at t=0)
[~, A_index(1)] = min(abs(grid.a - 1));
sim.A(1) = grid.a(A_index(1)) * (1 + shock_size);
sim.K(1) = ss.k;

% 4) Loop
for t = 1:T-1
    j  = A_index(t);
    At = sim.A(t);
    Kt = sim.K(t);

    Kt_eval = min(max(Kt, grid.k(1)), grid.k(end));  % clamp for evaluation

    % policies at (A_t,K_t)
    sim.C(t) = F_c(At, Kt_eval);
    sim.N(t) = F_n(At, Kt_eval);
    sim.I(t) = F_i(At, Kt_eval);
    sim.Y(t) = prod_function(At, Kt, sim.N(t), par);

    % capital via k'-policy
    sim.K(t+1) = F_k(At, Kt_eval);

    % --- stochastic Markov evolution for A ---
    p = P(j,:);                        % correct row
    u = rand;                          % shared RNG across runs controls GIRF
    cdf = cumsum(p);
    A_index(t+1) = find(u <= cdf, 1, 'first');
    sim.A(t+1)   = grid.a(A_index(t+1));
end

% 5) Last period
AT = sim.A(T); KT = sim.K(T);
KT_eval = min(max(KT, grid.k(1)), grid.k(end));
sim.C(T) = F_c(AT, KT_eval);
sim.N(T) = F_n(AT, KT_eval);
sim.I(T) = F_i(AT, KT_eval);
sim.Y(T) = prod_function(AT, KT, sim.N(T), par);

% 6) Moments (optional)
logY = log(sim.Y); logC = log(sim.C); logN = log(sim.N); logI = log(sim.I);
m.C=mean(sim.C); m.N=mean(sim.N); m.Y=mean(sim.Y); m.I=mean(sim.I);
sd.Y=std(logY); sd.C=std(logC); sd.N=std(logN); sd.I=std(logI);
logcorr.CY=corr(logC,logY); logcorr.NY=corr(logN,logY); logcorr.IY=corr(logI,logY);

% 7) Plots (optional)
fig1=figure('Color','w');
sgtitle(['Simulation (T=' num2str(T) '), \phi=' num2str(par.phi) ', shock=' num2str(shock_size)],'FontWeight','bold');
subplot(3,2,1); plot(sim.K(1:T_graph),'LineWidth',1.2); title('K'); xlabel t; ylabel level
subplot(3,2,2); plot(sim.C(1:T_graph),'LineWidth',1.2); title('C'); xlabel t; ylabel level
subplot(3,2,3); plot(sim.I(1:T_graph),'LineWidth',1.2); title('I'); xlabel t; ylabel level
subplot(3,2,4); plot(sim.N(1:T_graph),'LineWidth',1.2); title('N'); xlabel t; ylabel level
subplot(3,2,5); plot(sim.Y(1:T_graph),'LineWidth',1.2); title('Y'); xlabel t; ylabel level
set(gcf,'Position',[100 100 1000 600]);
if ~isempty(outdir), exportgraphics(fig1, fullfile(outdir, ...
   ['simulation_phi_' num2str(par.phi) '_shock_' num2str(shock_size) '.png']), 'Resolution',300); end
end
