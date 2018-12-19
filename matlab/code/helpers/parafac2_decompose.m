function [As, objfuncval] = parafac2_decompose(X, ncomps)
Options(3) = 0;
Options(5) = 1;
const = [0, 2];
%0 => no constraint,
%1 => orthogonality
%2 => nonnegativity
%3 => unimodality (and nonnegativitiy)
%4 => L1 fitting (will be imposed in all modes)
%5 => L1 fitting and nonnegativity (will be imposed in all modes)
[A,H,C,P, objfuncval]=parafac2(X, ncomps, const, Options);
As = cell(1, 4);
As{1} = A;
As{2} = H;
As{3} = C;
As{4} = P;

end