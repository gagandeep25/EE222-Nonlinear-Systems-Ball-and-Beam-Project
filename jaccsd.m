function [z, A] = jaccsd(fun, x)
% JACCSD Jacobian through complex step differentiation
% [z, A] = jaccsd(fun, x)
% z = fun(x)          : the function value at x
% A = fun'(x)         : the Jacobian matrix at x

z = fun(x);         % Evaluate the function at x to get z
n = numel(x);       % Get the number of elements in the input vector x
m = numel(z);       % Get the number of elements in the output vector z
A = zeros(m, n);    % Initialize the Jacobian matrix (m x n)

h = n * eps;        % Define the perturbation size (scaled by the size of x)

for k = 1:n
    x1 = x;                     % Copy the input vector x
    x1(k) = x1(k) + h * i;      % Perturb the k-th element of x by a small imaginary value
    A(:, k) = imag(fun(x1)) / h; % Compute the partial derivative by taking the imaginary part
end