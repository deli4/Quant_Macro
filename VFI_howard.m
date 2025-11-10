
function [v,c_pol,k_pol,n_pol,i_pol, grid,prob_mat]= VFI_howard (ss,par,N,M, tol, max_it,howard_steps)
% Value Function Iteration with the Howard method 
% We maximize over the bellman for the first 100 iterations to be on the
% right path, an dthen perform the policy iteration, while maximizing every
% howard_steps time. 
% ss : the steady state values 
% par : the parameters 
% N: the nb of points on the grid of K 
% M: nb of points on the grid of A 
% tol : tolerance level 
% max_it: maximum iterations of the VFI
% howard_steps: every X steps we perform the full maximization


%1. Define an equally spaced grid, 50 steps K1=0.9k* and K50=1.1k*
cpar.kmin=0.8*ss.k;
cpar.kmax=1.2*ss.k;

grid.k=linspace(cpar.kmin,cpar.kmax,N);

%2. Discreize the log(At) process - get a transition matrix and a grid on A
%Discretize the process for log(At) using Rouwenhorst's method (you can find the code at tobiasbroer.eu/teaching)
par.mu=0;
[grid.loga, prob_mat] =rouwenhorst(M,par.mu, par.rho, par.sigma); %log_A is the values of the discretized AR(1) process and prob_mat the transition matrix 


%%%%%%%%%%%%%%%%%
%vraiment unsure de ca !!!
grid.a=exp(grid.loga); %because our trasnition is in log 
%%%%%%%%%%%%%%%%%


% 3. Choose an initial value function v0(Ki,zj), an Na ×M matrix.
v_ini=zeros(M,N);

% initialize the policy function matrices and next guess functions 

c_pol=zeros(M,N);
k_pol=zeros(M,N);
i_pol=zeros(M,N);
n_pol=zeros(M,N);

% VFI 

iter=0;
err=inf;
v=v_ini;

tic 
while iter<max_it && err>tol
    iter=iter+1;
 
    v_next=zeros(M,N);

    if iter<100 || mod(iter,howard_steps)==1
      %%%%%%%%%%%%%%%%%%%%% MAXIMIZATION STEP %%%%%%%%%%%%%%%%%%%%%

       
        for j=1:M %where we loop over the different a
           a=grid.a(j);
           for i=1:N % loop over different K
                k=grid.k(i);
    
            
                %have the right hand side as a function c 
                bell_rhs= @(C) bellman_rhs(C,k,a,par,j,prob_mat,grid,v,M); 
    
                % maximize over c using the golden search
    
                %%%%%%%%%%%%%%%%
                % revoir le feasibility bound de c_upper 
                %%%%%%%%%%%%%%%%%
    
                c_lower = 1e-8;
                y_ub    = prod_function(a, k, 1, par);   % output at n=1
                c_upper = max(10*c_lower, 0.99*y_ub);
    
    
                % if ~(isfinite(c_upper) && c_upper > c_lower)
                %     v_next(i,j) = -1e12;  % kill state
                %     continue
                % end
    

    
                [c_max, v_next(j,i)] = golden_search_improved(bell_rhs, c_lower, c_upper, 1e-8, 1000); %chose the same tolerance and iterations as in VFI
    
    
                % store the results in the policy functions 
                c_pol(j,i)=c_max;
                n_pol(j,i)=labour_foc(c_max,k,a,par);
                y=prod_function(a,k,n_pol(j,i),par);
                i_pol(j,i)=y-c_max;
                k_pol(j,i) = lom_function(i_pol(j,i), k, par);
    
    
    
            end 
        end
    else 
        %%%%%%%%%%%%%%%%%%%%% POLICY FUNCTION ITERATION %%%%%%%%%%%%%%%%%%%%%
        for j=1:M
            for i=1:N
                c = c_pol(j,i);
                n = n_pol(j,i);
                k = grid.k(i);
                a = grid.a(j);


                y = prod_function(a, k, n, par);
                i_val = y - c;
                k_next = lom_function(i_val, k, par);


                %interpolate next step
                V_interp = zeros(M,1);
                for m = 1:M
                    V_interp(m) = linear_interpolation(grid.k, v(m,:), k_next);
                end
                EV = prob_mat(j,:) * V_interp;

                v_next(j,i) = utility_function(c, n, par) + par.beta * EV;
            end
        end
    end


% convergence 
err = max(abs(v_next(:) - v(:)));
v   = v_next;
            
if iter<50 || mod(iter,100)==1 
    disp(['iteration: ', num2str(iter), ' ; V(n+1)-V(n): ', num2str(err)])
end 


end 

toc


end
