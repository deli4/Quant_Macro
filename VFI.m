
function [v,c_pol,k_pol,n_pol,i_pol, grid, prob_mat]= VFI (ss,par,N,M, tol, max_it)

% Value Function Iteration
% ss : the steady state values 
% par : the parameters 
% N: the nb of points on the grid of K 
% M: nb of points on the grid of A 
% tol : tolerance level 
% max_it: maximum iterations of the VFI


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
v=zeros(M,N);
v_next=zeros(M,N);
      
% initialize the policy function matrices and next guess functions 

c_pol=zeros(M,N);
k_pol=zeros(M,N);
i_pol=zeros(M,N);
n_pol=zeros(M,N);

% VFI 

iter=0;
err=inf;


tic 
while iter<max_it && err>tol
    iter=iter+1;
 


     for j=1:M %where we loop over the different a
            a=grid.a(j);

            for i=1:N % loop over different K
                 k=grid.k(i);
    
                global penalties
                penalties = struct('C',0,'N',0,'Kneg',0,'Kout',0);
    
    
            
                %have the right hand side as a function c 
                bell_rhs= @(C) bellman_rhs(C,k,a,par,j,prob_mat,grid,v,M); 
    
                % maximize over c using the golden search
    
                c_lower = 1e-6; 
                y_ub    = prod_function(a, k, 1, par);   % output at n=1
                c_upper = y_ub-c_lower;
    
               
    
                [c_max, v_next(j,i)] = golden_search_improved(bell_rhs, c_lower, c_upper, 1e-8, 1000); %chose the same tolerance and iterations as in VFI
    
    
                % store the results in the policy functions 
                c_pol(j,i)=c_max;
                n_pol(j,i)=labour_foc(c_max,k,a,par);
                y=prod_function(a,k,n_pol(j,i),par);
                i_pol(j,i)=y-c_max;
                k_pol(j,i) = lom_function(i_pol(j,i), k, par);
    
    
                v_next(j,i)=v_next(j,i);
               



        end 
    end
     % display(['penalties on C: ',num2str(penalties.C), newline, ...
     %            'penalties on N: ', num2str(penalties.N), newline, ...
     %           'penalties on K prime because negative: ', num2str(penalties.Kneg), newline, ...
     %           'penalties on K prime because out of bounds: ', num2str(penalties.Kout)]);
% convergence 
err = max(abs(v_next(:) - v(:)));
v   = v_next;
            
if iter<50 || mod(iter,100)==1 
    disp(['iteration: ', num2str(iter), ' ; V(n+1)-V(n): ', num2str(err)])
end 

end 

toc


end
