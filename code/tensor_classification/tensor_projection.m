function Xprojections = tensor_projection(Xs, Us)

nsamples = length(Xs);
nmodes = length(Us);
origdims = size(Xs{1});
Xprojections = cell(nsamples, 1);

for isample = 1:nsamples
    Xprojection = Xs{isample};
    
    for jmode = 1:nmodes % do multiplication between Xi tensor and U matrices
        Xprojection = tmult(Xprojection, Us{jmode}', jmode);
    end
    Xprojections{isample} = Xprojection;
    
    
    if false % equivalent to using tmult, can be used if tmult is not available
        dims = origdims;
        for jmode = 1:nmodes % do multiplication between Xi tensor and U matrices
            Xprojection = Us{jmode}'*matricizing(Xprojection, jmode);
            dims(jmode) = size(Us{jmode}, 2);
            Xprojection = unmatricizing(Xprojection, jmode, dims);
            Xprojections{isample} = Xprojection;
        end;
    end
    
    
end

end