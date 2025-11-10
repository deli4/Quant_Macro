function bellman_rhs=bellman_rhs(C,K,A,par,j, transition_matrix, grid, V,M)
% Calculate the RHS of the Bellman equation
%
%Inputs : 
% C: consumption guess 
% K, A: come from the loop in the VFI these are the points in the grid 
% par: the parameters of calibration
% cpar: the parameters of the VFI (matrix dimension, etc.)
% j: the current shock Aj in the state-space 
% transition_matrix: a transition matrix found with Tauchen or Rouwenhorst
% grid: the grids on k
%V: the value matrix 
%
% Note: 
% we look at the feasibility constraints on N,C and K', if they aren't
% repsected we impose a penalty (bellman=0), and we use the following functions
% N=h(c,k)
%Y=f(k,n)
% I= Y-C


penalty=-1e10;
global penalties

% Check feasibility: penalty if C<=0
if C<=penalty
    bellman_rhs=penalty;
    penalties.C=penalties.C+1;
    return
end 

%compute N
N=labour_foc(C,K,A,par); 

%Check feasibility: penalty if N not in [0,1]
if N<=0 || N>1
    bellman_rhs=penalty;
     penalties.N=penalties.N+1;
    return
end

% Compute Y and I
Y=prod_function(A,K,N,par);

I=invest_fct(Y,C);


if I <= 0.1*Y
   bellman_rhs = penalty;

    return
end



%adding a lower bound on I made the bellman penalized and stop the
%iteration after on step


%Compute K'
K_prime=lom_function(I,K,par);

%check feasibility: penalty if K'<0 and that it stays on the grid 
if K_prime<=1e-8 
    bellman_rhs = penalty; 
    penalties.Kneg=penalties.Kneg+1;
    return
 % elseif  K_prime > 
 % grid.k(end) || K_prime < grid.k(1) 
 %     bellman_rhs = penalty; 
 %     penalties.Kout=penalties.Kout+1;
 %     return

end 




% caculate the interpolated value for V'(K')
%V_interp = zeros(1,M);
EV=0;

for m = 1:M
    %V_interp(m) = interp1(grid.k, V(:,m), K_prime, 'linear','extrap');
    % V_interp(m)=interp_kp(grid.k, V(:,m), K_prime);
    V_interp=linear_interpolation(grid.k, V(m,:), K_prime);
    EV=EV+transition_matrix(j,m)*V_interp;
end

% EV = transition_matrix(j,:) * V_interp.';


%calculate the utlity function
u=utility_function(C,N,par);

bellman_rhs=u+par.beta*EV;

end 