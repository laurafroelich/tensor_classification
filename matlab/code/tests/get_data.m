function [Xs, classes] = get_data(nsamples, xdims)

    if nargin < 2
        xdims = [5, 7];

        if nargin < 1
            nsamples = 500;
        end
    end

    classes = round(rand(nsamples,1))+1;
    Xs = cell(1, nsamples);

    for isample = 1:nsamples
        if classes(isample) == 1
            Xs{isample} =  randn(xdims);
        else

            Xs{isample} =  randn(xdims) + 10;
        end
    end 
end