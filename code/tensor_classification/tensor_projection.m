function Xprojections = tensor_projection(Xs, Us)

nmodes = length(Us);

is_cell_data = iscell(Xs);

if is_cell_data
    nsamples = length(Xs);
    Xprojections = cell(nsamples, 1);
else
    Xdim = size(Xs);
    nsamples = Xdim(end);
    y_mat_size = NaN(1, nmodes+1);
    for imode = 1:nmodes
        Udims = size(Us{imode});
        y_mat_size(imode) = Udims(1);
    end
    y_mat_size(nmodes+1) = nsamples;
    Xprojections = NaN(y_mat_size);
    Xs = permute(Xs, [nmodes+1, 1:nmodes]);
    Xprojections = permute(Xprojections, [nmodes+1, 1:nmodes]);
end

for isample = 1:nsamples
    if is_cell_data
        Xprojection = Xs{isample};
    else
        Xprojection = Xs(isample,:);
    end
    
    for jmode = 1:nmodes % do multiplication between Xi tensor and U matrices
        Xprojection = tmult(Xprojection, Us{jmode}', jmode);
    end
    if is_cell_data
        Xprojections{isample} = Xprojection;
    else
        Xprojections(isample,:)=Xprojection;
    end
    
    if ~is_cell_data
        Xprojections = permute(Xprojections, [2:(n_modes+1), 1]);
    end
    
end

end