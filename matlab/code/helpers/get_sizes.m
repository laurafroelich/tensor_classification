function [n_items, mode_sizes, n_modes] = get_sizes(Xs, varargin)
% [n_items, mode_sizes, n_modes] = get_sizes(Xs, varargin)
%
% Get number of items and sizes of modes.
%
% Input:
%
% Xs: Array with at least two modes.
% varargin: If given, may only contain one element, which must be an
%           integer specifying the index at which observations run along.
%           Default is 1.
%
% Output:
%
% n_items: Number of observations.
% mode_sizes: Length of each mode.
% n_modes: Number of modes.
if length(varargin) == 1
    observation_mode = varargin{1};
else
    observation_mode = 1;
end


all_mode_sizes = size(Xs);
all_modes = 1:length(all_mode_sizes);
non_observation_modes = all_modes;
non_observation_modes(observation_mode) = [];

mode_sizes = all_mode_sizes(non_observation_modes);
n_items = all_mode_sizes(observation_mode);
n_modes = length(mode_sizes);
end
