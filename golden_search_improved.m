function [x_max, f_max] = golden_search_improved(func, x_lower, x_upper, tol, max_iter)
% Golden-section search for a maximum on [x_lower, x_upper]

% enforce order
if x_lower > x_upper
    t = x_lower; x_lower = x_upper; x_upper = t;
end

r = 0.5*(sqrt(5)-1);    % golden ratio conjugate

xl = x_lower; 
xu = x_upper;

% interior points
d  = r*(xu - xl);
x1 = xu - d;
x2 = xl + d;

f1 = func(x1);
f2 = func(x2);

iter = 0;
while abs(xu - xl) > tol && iter < max_iter
    iter = iter + 1;

    if f1 > f2
        % maximum in [xl, x2]
        xu = x2;
        x2 = x1;   % reuse
        f2 = f1;
        d  = r*(xu - xl);
        x1 = xu - d;
        f1 = func(x1);  % one new eval
    elseif f1 < f2
        % maximum in [x1, xu]
        xl = x1;
        x1 = x2;   % reuse
        f1 = f2;
        d  = r*(xu - xl);
        x2 = xl + d;
        f2 = func(x2);  % one new eval
    else
        % equal values: shrink symmetrically without collapsing interval to a point
        xl = x1;
        xu = x2;
        d  = r*(xu - xl);
        x1 = xu - d;
        x2 = xl + d;
        f1 = func(x1);
        f2 = func(x2);
    end
end

% choose by function value, not position
if f1 >= f2
    x_max = x1; f_max = f1;
else
    x_max = x2; f_max = f2;
end
end
