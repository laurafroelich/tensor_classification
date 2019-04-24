function [cmean_m_xmeans, xi_m_cmeans, nis] = classbased_differences(Xs, classes)
% [cmean_m_xmean, xi_m_cmean, nis] = classbased_differences(Xs, classes)
%
% input:
% Xs: cell array of (multi-dimensional) matrices
% classes: classes of matrices
%
% output:
% cmean_m_xmean: class means minus overall mean
% xi_m_cmean: observations minus corresponding class mean

nsamples = length(Xs);
nclasses = length(unique(classes));

if isa(Xs, 'cell')
    Xs = cell_array_to_nd_array(Xs);
end

Xmean = mean(Xs, 1);

Xmeansclasses = cell(1, nclasses);
Xsumsclasses = cell(1, nclasses);
nis = NaN(1, nclasses);
for iclass = 1:nclasses
    inds = find(classes==iclass);
    Xsumsclasses{iclass} = Xs{inds(1)};
    for iind = 2:length(inds)
        Xsumsclasses{iclass} = Xsumsclasses{iclass} + Xs{inds(iind)};
    end
    nis(iclass) = length(inds);
    Xmeansclasses{iclass} = Xsumsclasses{iclass}/nis(iclass);
end

xi_m_cmeans = cell(1, nsamples);
for isample = 1:nsamples
    xi_m_cmeans{isample} = Xs{isample}-Xmeansclasses{classes(isample)};
end

cmean_m_xmeans = cell(1, nclasses);
for iclass = 1:nclasses
    cmean_m_xmeans{iclass} = Xmeansclasses{iclass}-Xmean;
end

end
