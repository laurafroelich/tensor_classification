function scaled_class_observations = classwise_scalar_multiply(class_observations, class_scalars)
% scaled_class_observations = classwise_scalar_multiply(class_observations, class_scalars)
%
% Multiply all scalar entries for each class by scalar. Observations must
% be placed along last mode.
%
% Mandatory input:
% class_observations: Array containing observations along last mode.
% class_scalars: Vector of length equal to the size of the last mode of
%                class_observations.
%
% Output:
% scaled_class_observations: Array of size equal to the input
% class_observations, but with entries with the same last index scaled by
% the same scalar, as specified in the input class_scalars.
[~, sizeobs, nmodes] = get_sizes(class_observations, length(size(class_observations)));
classwise_n_obs_sqrts = repmat(class_scalars, [sizeobs(1), 1, sizeobs(2:end)]);
permuted_classwise_n_obs_sqrts = permute(classwise_n_obs_sqrts, ...
    [1 3:(nmodes+1) 2]);

scaled_class_observations = class_observations.*permuted_classwise_n_obs_sqrts;

end