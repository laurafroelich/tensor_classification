function scaled_class_observations = classwise_scalar_multiply(class_observations, class_scalars)
[~, sizeobs, nmodes] = get_sizes(class_observations, length(size(class_observations)));
classwise_n_obs_sqrts = repmat(class_scalars, [sizeobs(1), 1, sizeobs(2:end)]);
permuted_classwise_n_obs_sqrts = permute(classwise_n_obs_sqrts, ...
    [1 3:(nmodes+1) 2]);

scaled_class_observations = class_observations.*permuted_classwise_n_obs_sqrts;

end