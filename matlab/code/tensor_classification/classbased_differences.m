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

n_modes = length(size(Xs)) - 1;

Xmean = mean(Xs, 1);

Xmeansclasses = cell(1, nclasses);
nis = NaN(1, nclasses);
% https://se.mathworks.com/matlabcentral/answers/367188-how-to-filter-an-n-d-array-with-matrix-indexing
C = repmat({':'},1,n_modes + 1); % colons for all trailing dimensions
xi_m_cmeans = Xs;

for iclass = 1:nclasses
    inds = find(classes==iclass);
    C{1} = inds; 
    
    nis(iclass) = length(inds);
    
    Xmeansclasses{iclass} = mean(subsref(Xs,substruct('()',C)), 1); 
    current_class_mean_obs_differences = ...
        subsref(Xs,substruct('()',C)) - Xmeansclasses{iclass};
    
    xi_m_cmeans = subsasgn(xi_m_cmeans, substruct('()',C),...
        current_class_mean_obs_differences); 
end

cmean_m_xmeans = cell(1, nclasses);
for iclass = 1:nclasses
    cmean_m_xmeans{iclass} = Xmeansclasses{iclass}-Xmean;
end

end
