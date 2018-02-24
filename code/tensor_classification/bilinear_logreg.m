function [f, g, H] = bilinear_logreg(w0uv, Xs, targets, D)
% Xs: 3D matrix
% w0uv: column vector

[I, J, ntrials] = size(Xs);
w0 = w0uv(1);
U = reshape(w0uv(2:(2+I*D-1)),I,D);
V = reshape(w0uv((2+I*D):end),J,D);
XsU = tmult(Xs, U', 1);
XsV = tmult(Xs, V', 2);
psiXssq = squeeze(sum(sum(bsxfun(@times, XsU, V'),1),2)); % psiXs squeezed

linfunc = psiXssq + w0; % vector of length ntrials
expectedlabels = 1./(1+exp(-linfunc));

loglikelihoods = targets.*linfunc - log(1+exp(linfunc));
f = sum(loglikelihoods);

trueminuspreds = targets - expectedlabels; % must be column vector

% Gradient:
dfdw = sum(trueminuspreds);
dfdu = tmult(XsV,trueminuspreds', 3);
dfdvt = tmult(XsU,trueminuspreds', 3)';

g = [dfdw; dfdu(:); dfdvt(:)];

f = -f;
g = -g;

if nargout > 2
% Hessian
H = NaN(1+(I+J)*D);
pioneminuspi = expectedlabels.*(1-expectedlabels);
d2fdw0dw0 = -sum(pioneminuspi);
d2fdw0dU = -tmult(XsV,pioneminuspi', 3);
d2fdw0dV = -tmult(permute(XsU,[2 1 3]),pioneminuspi', 3);
XsV3=matrizicing(XsV,3);
d2fdudu = -bsxfun(@times,XsV3,pioneminuspi)'*XsV3;
XsU3=matrizicing(permute(XsU,[2 1 3]),3);
d2fdvdvtemp = -bsxfun(@times,XsU3,pioneminuspi)'*XsU3;
d2fdvdv = d2fdvdvtemp;%reshape(d2fdvdvtemp, [D J ])
d2fdudv = kron(eye(D), tmult(Xs,trueminuspreds',3))'-bsxfun(@times,XsU3,pioneminuspi)'*XsV3;

H(:,1)=[d2fdw0dw0; d2fdw0dU(:); d2fdw0dV(:)];
H(2:end,2:(2+I*D-1))= [d2fdudu; d2fdudv]; 
H((2+I*D):end,(2+I*D):end)= d2fdvdv;

H=tril(H)+tril(H,-1)';

H = -H;
end

