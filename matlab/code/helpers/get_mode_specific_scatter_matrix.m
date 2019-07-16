function mode_specific_scatter_matrix = get_mode_specific_scatter_matrix(...
    data_difference_tensor, kmode, lowerdims, Us, n_obs, nmodes, sizeX)
% mode_specific_scatter_matrix = get_mode_specific_scatter_matrix(...
%    data_scatter_tensor, kmode, lowerdims, Us, n_obs, nmodes, sizeX)
%
% Mandatory input:
% data_difference_tensor: Tensor (multi-modal array) containing the 
%                         differences to base the scatter on.
% kmode: Integer, mode to get the scatter matrix for. The data differences 
%        will be projected unto all other modes before multiplying to get 
%        the squares of differences, i.e. scatters.
% lowerdims: Vector containing integers, with length equal to number of 
%            observation modes, i.e. length of Us, or, equivalently, 
%            length of size of the input data_difference_tensor minus one.
%            Each integer denotes the number of dimensions that
%            corresponding mode is projected to.
% Us: Cell array of projection matrices.
% n_obs: Integer, number of observations.
% nmodes: Integer, number of modes.
% sizeX: Vector containing integers of length equal to Us, size of each mode.
%
% Output:
% mode_specific_scatter_matrix: Matrix containing the projected squared
%                               differences with the mode specified in the
%                               input, kmode, along the first dimension.
permute_vector = 1:(nmodes+1);
permute_vector(1) = kmode;
permute_vector(kmode) = 1;

othermodes = setdiff(1:nmodes, kmode);

for othermode = othermodes
    data_difference_tensor=tmult(data_difference_tensor,Us{othermode}', othermode);
end

projected_data_differences_matrix=reshape(...
    permute(data_difference_tensor, permute_vector),...
    [sizeX(kmode), prod(lowerdims(othermodes))*n_obs]);
mode_specific_scatter_matrix = projected_data_differences_matrix*projected_data_differences_matrix';
end
