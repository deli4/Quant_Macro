function yi=linear_interpolation(x,y,xq)
% Linear Interpolation
%
% This function allows for extrapolation on the right and the left of the
% boundaried 
%
% INPUTS: 
% x: grid points- assorted in ascending order
% y: function value at grid points
% xq: query points where interpolation/extrapolation is needed 
%
% OUTPUT: y(xq), value of the function at xq



% check grid points in ascending 
if ~isvector(x)
    disp('Not a grid vector');
    return
end

% must be strictly increasing
if any(diff(x) <= 0)
    disp('not sorted or not a grid vector');
    return
end

% Check xq is a scalar 
if isscalar(xq)
    xq_vec=xq;
else 
    xq_vec=xq(:);
end 

n=length(xq_vec);
yi=zeros(size(xq_vec));

for k=1:n
    xq=xq_vec(k);

   
    if xq<=x(1)  %1. extrapolation on the left
        if xq>0 %k' is out of the grid but we need to make sure it respects feasibilty constraint k'>0
            slope=(y(2)-y(1))/(x(2)-x(1));
            yi(k)=y(2)+slope*(xq-x(2));

            %this additional penalty makes the VFI diverge 
        else % if k' does not respect feasibility constraint
            yi(k)=-1e12; %penalty
            printf('Warning: Left extrapolation at xq = %.4f --> Penalize\n', xq);
        end 
   
    elseif xq>=x(end) %2. extrapolation on the right
        slope=(y(end)-y(end-1))/(x(end)-x(end-1));
        yi(k)=y(end)+slope*(xq-x(end));
    
    else % 3.interpolation
      % Interpolation: find the interval [x(idx), x(idx+1)]
    idx = find(x <= xq, 1, 'last');
    
    % Linear interpolation formula
    weight = (xq - x(idx)) / (x(idx+1) - x(idx));
    yi(k) = y(idx) + weight * (y(idx+1) - y(idx));
    end 
end 

% Reshape output to match input shape
    if ~isscalar(xq)
        yi = reshape(yi, size(xq));

end 



