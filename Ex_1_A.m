%% PS 2 - Exercise 1 - Value Function Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Délia Manso de Zuniga 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Establish all paths 
%set directory of the main code
%function path 
function_path= addpath(pwd,'01_functions');

% graph path 
graph_save = fullfile(pwd,'02_graphs'); %the graph folder


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Calibration and VFI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. find the value of theta 

% declare the parameters 
par.alpha=0.36;
par.beta=0.99;
par.delta=0.025;
par.rho=0.95;
par.sigma=0.007;
par.gamma=1.0000001;
par.varphi=1;
par.phi=5;

% K/Y= beta*alpha/(1-beta*(1-beta))
ss.k_y=par.beta*par.alpha/(1-par.beta*(1-par.delta));
%N/Y=(K/Y)^(-alpha/(1-alpha))
ss.n_y=ss.k_y^(-par.alpha/(1-par.alpha));

%when n=1/3 we have theta: 
n=1/3;
par.theta=n^(-(par.varphi+par.gamma))*((1/ss.n_y)/(1/ss.k_y))^(-par.gamma)*(1-par.alpha)*((1/ss.k_y)-par.delta)^(-par.gamma)*(1/ss.n_y);


% the steady state value of output is 
ss.y=ss.k_y^(par.alpha/(1-par.alpha))*n;

%the steady state value of capital is
ss.k=ss.y*ss.k_y;

%the steady state value of consumption is 
ss.c=ss.y-par.delta*ss.k;

disp(['the steady state value of the ratio of capital to output, K/Y=',num2str(ss.k_y),newline, ...
    'the steady state value of the ratio of labour to output, N/Y=', num2str(ss.n_y),newline,  ...
    'the steady state value of theta when N=1/3, theta= ', num2str(par.theta), newline, ...
    'the steady state value of output when N=1/3, Y*=', num2str(ss.y), newline, ...
    'the steady state value of K when N=1/3, K*=', num2str(ss.k), newline ,...
    'the steady state value of C when N=1/3, C*=', num2str(ss.c)
    ])


%% 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  VFI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=50; % 50 steps
M=5;
tol=1e-6;
max_it=2000;


% phi = 0
par.phi = 0;
[v0, c0, k0, n0, i0, grid, P] = VFI(ss, par, N, M, tol, max_it);

% phi = 5
 par.phi = 5;
[v5, c5, k5, n5, i5, grid, P] = VFI(ss, par, N, M, tol, max_it);

% plot
plot_VFI_results(grid, par, v0, v5, c0, c5, k0, k5, n0, n5, i0, i5, graph_save)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Howard improvement algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It takes 40 sec to run for both values of phi and yields the same result,
% so for efficiency, we run this algorithm
% 
% par.phi = 0;
% howard_steps=50; 
% [v0_h, c0_h, k0_h, n0_h, i0_h, grid,P] = VFI_howard(ss, par, N, M, tol, max_it,howard_steps);
% 
% par.phi = 5;
% [v5_h, c5_h, k5_h, n5_h, i5_h, grid,P] = VFI_howard(ss, par, N, M, tol, max_it,howard_steps);
% 
% % PLOT AND SAVE GRAPH in .dir/figures
% plot_VFI_results(grid, par, v0_h,v5_h , c0_h, c5_h, k0_h, k5_h, n0_h, n5_h, i0_h, i5_h,graph_save)

 %% VFI with a spline 
% 
% par.phi = 0;
% [v0s, c0s, k0s, n0s, i0s, grid, P] = VFI_spline(ss, par, N, M, tol, max_it);
% 
% opts0s = struct('drop_first_row',true, 'title_suffix', 'spline, \phi = 0', 'cols', [1 M]);
% h_zero  = plot_policies_by_A(grid.k(:), grid.a(:), c0s, n0s, i0s, k0s, opts0s);
% hV_zero = plot_value_by_A  (grid.k(:), grid.a(:), v0s, opts0);

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.phi=0;
marg_u_c= @(C) C^(-par.gamma); %marginal utility of consumption  

accuracy(par, grid, c0, n0,k0, P, marg_u_c);

par.phi=5;

accuracy(par, grid, c5, n5, k5,P, marg_u_c);

%print out the euler errors statistics 
% mean, sd 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1000 iterations and random shocks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% baseline path for the IRF- at t=1, economy is at steady state 
no_shock=0;
H = 200; % plotting horizon 



%phi=0
par.phi=0;

[sim0,mean0,sd0,logcorr0]=simulation_irf(par,grid,P,ss,c0,n0,i0,k0,10000,12345, no_shock,H,graph_save);

disp([ 'FOR phi=0, the mean of', newline,...
    'log(C): ',num2str(mean0.C), newline, ...
    'log(N): ',num2str(mean0.N), newline, ...
    'log(I): ', num2str(mean0.I), newline, ...
    'log(Y): ' num2str(mean0.Y),newline]);

disp(['FOR phi = 0, the standard deviation of:', newline, ...
    'log(C): ', num2str(sd0.C), newline, ...
    'log(N): ', num2str(sd0.N), newline, ...
    'log(I): ', num2str(sd0.I), newline, ...
    'log(Y): ', num2str(sd0.Y),newline]);


disp(['FOR phi=0, the correlation to output is: ', newline, ...
    'corr(C,Y)= ', num2str(logcorr0.CY), ...
    'corr(N,Y)= ', num2str(logcorr0.NY), newline, ...
    'corr(I,Y)= ', num2str(logcorr0.IY),newline])

% compute the relative volatility
rel0.C=sd0.C/sd0.Y;
rel0.N=sd0.N/sd0.Y;
rel0.I=sd0.I/sd0.Y;

display(['FOR phi=0, the volatility of each input compared to output:', newline, ...
    'rel. vol C: ',num2str(rel0.C), newline, ...
    'rel vol N: ',num2str(rel0.N), newline, ...
    'rel vol I: ', num2str(rel0.I),newline])

%phi=5
par.phi=5;


[sim5,mean5,sd5,logcorr5]=simulation_irf(par,grid,P,ss,c5,n5,i5,k5,10000,12345,no_shock,H, graph_save);


disp([ 'FOR phi=5, the mean of', newline,...
    'log(C): ',num2str(mean5.C), newline, ...
    'log(N): ',num2str(mean5.N), newline, ...
    'log(I): ', num2str(mean5.I), newline, ...
    'log(Y): ' num2str(mean5.Y), newline]);


disp(['FOR phi = 5, the standard deviation of:', newline, ...
    'log(C): ', num2str(sd5.C), newline, ...
    'log(N): ', num2str(sd5.N), newline, ...
    'log(I): ', num2str(sd5.I), newline, ...
    'log(Y): ', num2str(sd5.Y), newline]);

disp(['FOR phi=5, the correlation to output is: ', newline, ...
    'corr(C,Y)= ', num2str(logcorr5.CY), ...
    'corr(N,Y)= ', num2str(logcorr5.NY), newline, ...
    'corr(I,Y)= ', num2str(logcorr5.IY), newline])

% compute the relative volatility
rel5.C=sd5.C/sd5.Y;
rel5.N=sd5.N/sd5.Y;
rel5.I=sd5.I/sd5.Y;

display(['FOR phi=5, the volatility of each input compared to output:', newline, ...
    'rel. vol C: ',num2str(rel5.C), newline, ...
    'rel vol N: ',num2str(rel5.N), newline, ...
    'rel vol I: ', num2str(rel5.I), newline])


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocked Simulation + IRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shock=0.01;

%phi=0
par.phi=0;

[shock0,mean_shock0,sd_shock0,logcorr_shock0]=simulation_irf(par,grid,P,ss,c0,n0,i0,k0,10000,12345, shock,H,graph_save);

%phi=5
par.phi=5;

[shock5,mean_shock5,sd_shock5,logcorr_shock5]=simulation_irf(par,grid,P,ss,c5,n5,i5,k5,10000,12345, shock,H,graph_save);


% --- IRF plots in levels ---
t = 0:H-1;
vars   = {'A','C','N','I','Y','K'};
titles = {'TFP','Consumption','Labor','Investment','Output','Capital'};

figIRF = figure('Color','w');     % <-- capture handle here
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

for k = 1:numel(vars)
    nexttile;
    % φ = 0 (blue) and φ = 5 (orange)
    plot(t, shock0.(vars{k})(1:H), 'b', 'LineWidth', 1.5); hold on;
    plot(t, shock5.(vars{k})(1:H), 'r', 'LineWidth', 1.5);
    % steady-state reference line (for all but A)
    if ~strcmp(vars{k}, 'A')
        yline(sim0.(vars{k})(1),'k:');
    else
        yline(1,'k:');  % steady state A=1
    end
    title([titles{k} ' (levels)']);
    xlabel('Time (periods)'); ylabel('Level');
    legend('\phi=0','\phi=5','Location','best'); legend boxoff;
    hold off;
end

% --- Global title and subtitle ---
sgtitle({'\bfIRF for a 1% increase in A at t=0', ...
          'for \phi = 0 and \phi = 5'}, 'FontSize', 12);

% --- Resize for readability ---
set(figIRF, 'Position', [100 100 1000 700]);   % use fig handle, not gcf

% --- Ensure folder exists ---
if ~exist(graph_save, 'dir')
    mkdir(graph_save);
end

% --- Save the figure ---
title_graph = 'IRF_phi0_phi5.png';
exportgraphics(figIRF, fullfile(graph_save, title_graph), 'Resolution', 300);

