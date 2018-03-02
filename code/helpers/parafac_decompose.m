function [As, objfuncval] = parafac_decompose(X, ncomps)
Options(2) = 2;  % Initialise with random values.
Options(5) = NaN; % how often to show fit
const = [0, 0, 1]; % Constraints for parafac. For tucker, the default of
% orthogonality is used.
%0 => no constraint,
%1 => orthogonality
%2 => nonnegativity
%3 => unimodality (and nonnegativitiy)
%4 => L1 fitting (will be imposed in all modes)
%5 => L1 fitting and nonnegativity (will be imposed in all modes)

[As, ~, objfuncval]=parafac(X, ncomps, Options, const);
end