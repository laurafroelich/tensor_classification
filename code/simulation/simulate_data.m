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

% the patterns should be the same for both classes
U1 = orth(randn([nrows, ncomps])); % simulate orthogonal matrix for mode 1
U2 = orth(randn([ncols, ncomps])); % simulate orthogonal matrix for mode 2

class1size = floor(nobs/2);
class2size = nobs - class1size;

% only the cores are assumed to differ systematically between the classes
%cores1 = repmat(randn([ncomps, ncomps])*sqrt(0.01), [1, 1, class1size]); %simulate_core_matrices(class1size, ncomps, parafac_structure);
%cores2 = repmat(randn([ncomps, ncomps])*sqrt(0.01), [1, 1, class2size]); %simulate_core_matrices(class2size, ncomps, parafac_structure);
cores1 = repmat(diag(diag(randn([ncomps, ncomps])))*sqrt(5), [1, 1, class1size]); %simulate_core_matrices(class1size, ncomps, parafac_structure);
cores2 = repmat(diag(diag(randn([ncomps, ncomps])))*sqrt(5), [1, 1, class2size]); %simulate_core_matrices(class2size, ncomps, parafac_structure);

% concatenate the cores for the two classes
cores = cat(3, cores1, cores2);
cores = cores + randn(size(cores))*sqrt(1); %sqrt(2) in latest sent plot

% generate the class labels for each observation
y1 = zeros(class1size, 1);
y2 = ones(class2size, 1);
y = cat(1, y1, y2)+1;

% shuffle the labels and the observations identically
shuffled_order = randperm(nobs);
y = y(shuffled_order);
cores = cores(:,:,shuffled_order); %+ randn(size(cores))*sqrt(0.1);
sigma  = 5;
Nnoisecomp = ncomps;
Unoise1 = orth(randn([nrows, Nnoisecomp])*sqrt(0.01)); % simulate orthogonal matrix for mode 1
Unoise2 = orth(randn([ncols, Nnoisecomp])*sqrt(0.01)); % simulate orthogonal matrix for mode 2

% generate the simulated observations by multiplying the cores and the
% factors for each mode, followed by noise addition.
x = tmult(...
    tmult(cores, U1, 1),...
    U2, 2);
x = x + sigma*tmult(tmult(randn(Nnoisecomp,Nnoisecomp,nobs), Unoise1, 1), Unoise2, 2);% + randn(size(x))*sqrt(0.01);% + ... % sqrt(0.1) is good. With 1, CMDA fails to estimate with NaN/Inf values.
    %... % with sqrt(0.8), sqrt(0.4) and sqrt(0.5), DATEReig fails. (when drawing cores form Wishart).
    %reshape(repmat(randn(1, size(x, 3))*sqrt(10), [nrows*ncols, 1]), [nrows, ncols, size(x,3)]); %randn(1)*sqrt(0);

end



function Z = posdef(n)
% Z = posdef(n)
% Generate a symmetric positive definite matrix of dimensions n by n.
% https://stat.ethz.ch/pipermail/r-help/2008-February/153708

% draw a vector of dimensions n by one from the standard uniform
% distribution on the interval 0 to 1.
eigen_vals = rand([n, 1]);

% draw n*n standard normally distributed numbers arranged in an n by n
% matrix.
Z = randn([n, n]);
[Q, R] = qr(Z);
d = diag(R);
ph = d ./ abs(d);
Omat = Q * diag(ph);
Z = Omat' * diag(eigen_vals) * Omat;
end

function cores = simulate_core_matrices(nobs, ncomps, parafac_structure)
% correlation matrix of the simulated cores, must be positive definite and
% symmetric to be a correlation matrix.
Z = posdef(ncomps);

% calculate R for speed-up when calling wishrnd multiple times.
R = chol(Z);

% lower numbers for the degrees of freedom imply higher variance for the
% Wishart distribution.
df = 1;

% draw nobs matrices from the Wishart distribution with correlation Z.
cores = NaN([ncomps, ncomps, nobs]);
for iobs = 1:nobs
    if ~parafac_structure
        cores(:, :, iobs) = wishrnd(Z, df, R);
    else % when the PARAFAC structure is simulated, the core matrix is 
        % diagonal.
        cores(:, :, iobs) = diag(diag(wishrnd(Z, df, R)));
    end
end


end
