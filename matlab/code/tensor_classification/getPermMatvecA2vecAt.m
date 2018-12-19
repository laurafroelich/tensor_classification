function Tmn = getPermMatvecA2vecAt(m, n)
% Tmn = getPermMatvecA2vecAt(m, n)
%
% Returns the permutation matrix Tmn such that the following holds for
% the m*n matrix A and n*m matrix At = A': Tmn*A(:) = At(:)

Tmn = zeros(m*n, m*n);
for in = 1:n
Tmn(((in-1)*n*m^2+in):(m*n+n):(m^2*in*n)) = 1;
end

% for testing
% A = randn(5,9);
%At = A';
%Tmn*A(:) == At(:)

end
