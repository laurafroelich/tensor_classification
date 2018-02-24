function [f, g] = bilinear_logreg_tucker(w0UVW, Xs, targets, Ks)
% Xs: 3D matrix
% w0uv: column vector
% Ks: vector of components for each mode

if length(Ks) == 1
   Ks = repmat(Ks, 1, length(size(Xs))-1); 
end
   K = Ks(1); % for simplicity, we only allow same number of components in each mode
   Ks(2) = K; % for simplicity, we only allow same number of components in each mode

[I, J, ~] = size(Xs);
w0 = w0UVW(1);
U = reshape(w0UVW(2:(2+I*Ks(1)-1)),I,Ks(1));
V = reshape(w0UVW((2+I*Ks(1)):(2+I*Ks(1)+J*Ks(2)-1)),J,Ks(2));
W=nan(K);
ignore=logical(eye(K));
W(~ignore) = w0UVW((2+I*Ks(1)+J*Ks(2)):end);
W(ignore)=1;
XsU = tmult(Xs, U', 1);
XsV = tmult(Xs, V', 2);
%psiXs = squeeze(sum(sum(bsxfun(@times, XsU, V'),1),2)); 
psiXs = squeeze(sum(sum(bsxfun(@times, tmult(XsU, V', 2), W),1),2)); 

linfunc = psiXs + w0; % vector of length ntrials
expectedlabels = 1./(1+exp(-linfunc));

loglikelihoods = targets.*linfunc - log(1+exp(linfunc));
f = sum(loglikelihoods);

trueminuspreds = targets - expectedlabels; % must be column vector

% Gradient:
dfdw0 = sum(trueminuspreds);
dfdU = tmult(XsV,trueminuspreds', 3)*W';
dfdVt = tmult(XsU,trueminuspreds', 3)'*W;
dfdW = tmult(tmult(XsU, V', 2),trueminuspreds', 3);
dfdWtmp=dfdW(~ignore);
g = [dfdw0; dfdU(:); dfdVt(:); dfdWtmp(:)];
f = -f;
g = -g;


