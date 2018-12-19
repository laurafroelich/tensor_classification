function [As, objfuncval] = tucker2_decompose(X, ncomps)
Options(2) = 2;  % Initialise with random values.
Options(5) = NaN; % how often to show fit
[As, ~, objfuncval]=tucker(X, [ncomps, ncomps, -1], Options);
end