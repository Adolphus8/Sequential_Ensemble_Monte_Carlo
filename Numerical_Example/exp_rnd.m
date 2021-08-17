function [output] = exp_rnd(mu,bounds,N,M)
%% Exponential distribution random number generator:
%
% This function-handle generates a N-by-M samples from an exponential
% distribution within a fixed bound.
%
%--------------------------------------------------------------------------
% Author:
% Adolphus Lye         - adolphus.lye@liverpool.ac.uk
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%
% Inputs:
% mu:     Scalar value of the mean parameter of the exponential distribution;
% bounds: A 1 x 2 vector of bounds within which to generate samples;
% N:      Scalar value of row entries;
% M:      Scalar value of column entries;
% 
% Output(s):
% output: N x M matrix of samples;
%
%--------------------------------------------------------------------------
%% Defining the function:

% Define the Box function:
bounds = sort(bounds);
box = @(x) unifpdf(x, bounds(1),bounds(2));

output = zeros(N,M);
for i = 1:M
for j = 1:N
    
    
while true
output(j,i) = exprnd(mu, 1);
samp = output(j,i);
if box(samp)
% The box function is the Uniform PDF in the feasible region.
% If a point is out of bounds, this function will return 0 = false.
break;
end
end

end
end

output = sort(output);
end

