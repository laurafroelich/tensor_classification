function diagonal_indices = get_diagonal_indices(k, varargin)
% diagonal_indices = get_diagonal_indices(k)
%
% Get linear indices of elements on multi-dimensional diagonal.
%
% Input: 
%      k : integer, dimensionality of multidimensional cube.
%
% varargin{1}: number of modes in cube. Default: 2.
%
% Output:
%      vector containing linear indices of diagonal

if isempty(varargin)
    n_modes = 2;
else
    n_modes = varargin{1};
end

diagonal_indices = zeros(k, 1);

mat = reshape(1:(k^n_modes), repelem(k, n_modes));

S.type = '()';
    
for ik = 1:k
    S.subs = repmat({ik}, 1, n_modes);
    diagonal_indices(ik, 1) = subsref(mat, S);
end

% attempts at vectorising function to improve speed
if false
if n_modes == k
    diagonal_indices = diag(cumsum([1:(k+1):k^2; k^2.*ones(k-1,k)]));
else
    if n_modes == 2
        diagonal_indices = (1:(k+1):k^2)';
    else
        diagonal_indices = diag(cumsum([1:(k^(n_modes-1)+1):k^n_modes; k*(k+1).*ones(k-1,k)]));
        
        % k: 2, n_modes: 4                                             3*2
        %diagonal_indices = diag(cumsum([1:(k^(n_modes-1)+1):k^n_modes; factorial(n_modes-(n_modes-k)+1).*ones(k-1,k)]));
        
        % k: 3, n_modes:4                                             3*2*2
        %diagonal_indices = diag(cumsum([1:(k^(n_modes-1)+1):k^n_modes; ((k*n_modes)).*ones(k-1,k)]));
        

    end
end

end
end