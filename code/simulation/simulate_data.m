function [x, y] = simulate_data(nobs, nrows, ncols, ncomps, varargin)
% [x, y] = simulate_data(nobs, nrows, ncols)
%
% Simulate data based on either a Tucker or a PARAFAC structure.
% For a Tucker model, each observation is simulated as
% x_i =  U1 * core * U2',
% where U1 is a matrix of size nrows by ncomps and U2 is a matrix of size
% ncols by ncomps. U1 and U2 are simulated orthogonal matrices, and are
% the same for all observations. The matrix core has size ncomps by ncomps
% and is simulated independently for each observation, drawn from the
% distribution for the observation's class. If the first element of
% varargin is true, then the core matrix is diagonal.

if isempty(varargin)
    parafac_structure = false;
else
    parafac_structure = varargin{1};
end

U1 = orth(randn([nrows, ncomps])); % simulate orthogonal matrix for mode 1
U2 = orth(randn([ncols, ncomps])); % simulate orthogonal matrix for mode 2

class1size = floor(nobs/2);
class2size = nobs - class1size;

cores1 = simulate_core_matrices(class1size, ncomps, parafac_structure);
cores2 = simulate_core_matrices(class2size, ncomps, parafac_structure);

cores = cat(3, cores1, cores2);

y1 = zeros(class1size, 1);
y2 = ones(class2size, 1);
y = cat(1, y1, y2)+1;

shuffled_order = randperm(nobs);
y = y(shuffled_order);
cores = cores(:,:,shuffled_order);
x = tmult(...
    tmult(cores, U1, 1),...
    U2, 2);



end



function Z = posdef(n)
% https://stat.ethz.ch/pipermail/r-help/2008-February/153708
eigen_vals = rand([n, 1]);
Z = randn([n, n]);
[Q, R] = qr(Z);
d = diag(R);
ph = d ./ abs(d);
Omat = Q * diag(ph);
Z = Omat' * diag(eigen_vals) * Omat;
end

function cores = simulate_core_matrices(nobs, ncomps, parafac_structure)
% correlation matrix of the simulated cores
Z = posdef(ncomps);
R = chol(Z);
df = 100;
cores = NaN([ncomps, ncomps, nobs]);
for iobs = 1:nobs
    if ~parafac_structure
        cores(:, :, iobs) = wishrnd(Z, df, R);
    else
        cores(:, :, iobs) = diag(diag(wishrnd(Z, df, R)));
    end
end


end
