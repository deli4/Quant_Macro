function [sim,m,sd,logcorr]=simulation(par,grid,P,ss,c_pol,n_pol,i_pol,k_pol,T,seed,shock_size,T_graph, outdir)
% Simulate path of variables over time
% we simulate a N(O,1) for a sequence of A(t), t=1,...,T
%
% INPUTS:
% par: calibration parameters (structure)
% grid: for a and k (structure)
% P: transition matrix of A
% ss: steady state values (structure)
% each policy functions : for C,N,I,K
% T: nb of periods 
% seed: for the random generation (we chose 12345),
% shock_size: size of the shock on A at t=1
% outdir: folder where to store the plot
%
%Output: 
% The simulated path for A,C,N,I,Y and K, U the uniforms used for A, A_index  (structure)
% mean, sd, logcorr: series statistics and correlation to output (in a
% structure)
% rel: related volatility of C,N,I compared to output 





% the random sequence 
rng(seed);

%1. Initialize the vectors 
sim.A=zeros(T,1);
sim.C=zeros(T,1);
sim.N=zeros(T,1);
sim.I=zeros(T,1);
sim.Y=zeros(T,1);
sim.K=zeros(T,1);
A_index=zeros(T,1); %store index from grid.a of the current state At
%A_index tells me in which row of the transition matrix I am in


% Matrices to store the first step of interpolation 
C_values0 = zeros(1, 50); 
I_values0 = zeros(1, 50); 
N_values0 = zeros(1, 50);

%build extrapolant outside the loop
F_c = griddedInterpolant({grid.a, grid.k}, c_pol, 'linear','nearest'); 
F_n = griddedInterpolant({grid.a, grid.k}, n_pol, 'linear','nearest');
F_i = griddedInterpolant({grid.a, grid.k}, i_pol, 'linear','nearest');
F_k = griddedInterpolant({grid.a, grid.k}, k_pol, 'linear','nearest');

%  init (steady state, one-time level shock)
%  before the loop (one-time shock at t=0 if needed) 
[~, A_index(1)] = min(abs(grid.a - 1));
sim.A(1) = grid.a(A_index(1)) * (1 + shock_size);   % impact only at t=0
sim.K(1) = ss.k;


for t = 1:T-1
    At = sim.A(t);
    j  = A_index(t);   
    Kt = sim.K(t);
    Kt_eval = min(max(Kt, grid.k(1)), grid.k(end)); %leave K on the grid

    % policies at (A_t, K_t)
    sim.C(t) = F_c(At, Kt_eval);
    sim.N(t) = F_n(At, Kt_eval);
    sim.I(t) = F_i(At, Kt_eval);
    % use k'-policy directly to update capital
    sim.K(t+1) = F_k(At, Kt_eval);

    % output
    sim.Y(t) = prod_function(At, Kt, sim.N(t), par);

    if shock_size>0

     % deterministic A-path= for IRF
     sim.A(t+1) = exp(par.rho * log(At));
    else
        % Uncomment for stochastic= random
    p = P(j, :);                               % row of transition matrix
    A_index(t+1) = randsample(numel(grid.a), 1, true, p);
    sim.A(t+1)   = grid.a(A_index(t+1));       % set next-period A ON the grid
    end 
end

% last period T
AT = sim.A(T); KT = sim.K(T);
KT_eval = min(max(KT, grid.k(1)), grid.k(end));
sim.C(T) = F_c(AT, KT_eval);
sim.N(T) = F_n(AT, KT_eval);
sim.I(T) = F_i(AT, KT_eval);
sim.Y(T) = prod_function(AT, KT, sim.N(T), par);


%% Compute moments
% take logs 
logY = log(sim.Y);
logC = log(sim.C);
logN = log(sim.N);
logI = log(sim.I);

%mean 
m.C=mean(sim.C);
m.N=mean(sim.N);
m.Y=mean(sim.Y);
m.I=mean(sim.I);

% st. deviation
sd.Y = std(logY);
sd.C = std(logC);
sd.N = std(logN);
sd.I = std(logI);


% correlation to output
logcorr.CY = corr(logC, logY);
logcorr.NY = corr(logN, logY);
logcorr.IY = corr(logI, logY);



%% graph the simulated paths 

% === Plot simulated time series ===
fig1=figure('Color','w');  % white (light) background
    sgtitle(['\bfSimulation for T = ' num2str(T) ' (\phi = ' num2str(par.phi) '- shock size ' num2str(shock_size) ')' ], ...
            'Color','k', 'FontWeight','bold', 'FontSize',12);
    
    subplot(3,2,1)
    plot(sim.K(1:T_graph), 'b', 'LineWidth',1.2)
    title('\bfCapital (K_t)', 'Color','k')
    xlabel('Time','Color','k')
    ylabel('Level','Color','k')
    set(gca,'Color','w')
    
    subplot(3,2,2)
    plot(sim.C(1:T_graph), 'b', 'LineWidth',1.2)
    title('\bfConsumption (C_t)', 'Color','k')
    xlabel('Time','Color','k')
    ylabel('Level','Color','k')
    set(gca,'Color','w')
    
    subplot(3,2,3)
    plot(sim.I(1:T_graph), 'b', 'LineWidth',1.2)
    title('\bfInvestment (I_t)', 'Color','k')
    xlabel('Time','Color','k')
    ylabel('Level','Color','k')
    set(gca,'Color','w')
    
    subplot(3,2,4)
    plot(sim.N(1:T_graph), 'b', 'LineWidth',1.2)
    title('\bfLabor (N_t)', 'Color','k')
    xlabel('Time','Color','k')
    ylabel('Level','Color','k')
    set(gca,'Color','w')
    
    subplot(3,2,5)
    plot(sim.Y(1:T_graph), 'b', 'LineWidth',1.2)
    title('\bfOutput (Y_t)', 'Color','k')
    xlabel('Time','Color','k')
    ylabel('Level','Color','k')
    set(gca,'Color','w')
    
    
    set(gcf, 'Position', [100 100 1000 600]);  % resize window
    title_graph=['simulation_phi_' num2str(par.phi) ' ,shock size ' num2str(shock_size) '.png'];
    exportgraphics(fig1, fullfile(outdir,title_graph), 'Resolution',300);

end 