function [Tmax, lowerdims, usestoppingcrit, Us] = parse_varargin(...
    sizeX, nmodes, varargin)


if length(varargin) >= 1 && ~isempty(varargin{1})
    Tmax = varargin{1};
else
    Tmax = 100;
end

if length(varargin)>=2 && ~isempty(varargin{2})
    lowerdims = varargin{2};
else
    lowerdims = sizeX;
end

if length(varargin)>=3 && ~isempty(varargin{3})
    usestoppingcrit = varargin{3};
else
    usestoppingcrit = true;
end

if length(varargin)>=4 && ~isempty(varargin{4})
    Us = varargin{4};
    if ischar(Us)
        switch Us
            case 'randinit'
                Us = cell(1, nmodes);
                for kmode = 1:nmodes
                    Us{kmode} = orth(randn(sizeX(kmode), lowerdims(kmode)));
                end
            case 'ones'
                Us = cell(1, nmodes);
                for kmode = 1:nmodes
                    Us{kmode} = eye(sizeX(kmode), lowerdims(kmode));
                end
            case 'identity'
                Us = cell(1, nmodes);
                for kmode = 1:nmodes
                    Us{kmode} = eye(sizeX(kmode), lowerdims(kmode));
                end
            otherwise
                warning('Initialisation method not recognised.')
        end
    end
end


end