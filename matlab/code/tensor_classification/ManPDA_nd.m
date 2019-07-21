function [Us, outputs, Ys] = ManPDA_nd(Xs, classes, varargin)
% [Us, outsmode1, outsmode2, outerwhile, Ys] = ManPDA(Xs, classes, varargin)
% classes: vector containing class labels. Classes must be sequential
% numbers starting from one.
% varargin{1}: lowerdims
% varargin{2}: Us
% varargin{3}: opts
% varargin{4}: usestoppingcrit
% varargin{5}: maxits

%rng(0)
% mymat1 = reshape(1:12, 2, 3, 2)
% mymat2 = reshape(13:24, 2, 3, 2)
% Xs{1} = mymat1
% Xs{2} = mymat2

Xsample1 = Xs{1};
sizeX = size(Xsample1);

nmodes = length(sizeX);


if nargout > 2
    [Us, outputs, Ys] = ManPDA(Xs, classes, varargin);
end

if nargout > 1
    [Us, outputs] = ManPDA(Xs, classes, varargin);
end


if nargout == 1
    for imode = 1:nmodes
        Us = ManPDA(Xs, classes, varargin);
    end
end

end

