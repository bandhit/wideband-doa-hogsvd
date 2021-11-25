%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Frequency Selection
%
% Description : Frequency Selection
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 12 November 2016, Bandhit Suksiri,
%               Updated: 15 January  2019, Bandhit Suksiri.
%
% Copyright 2016 - 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [head_idx, tail_idx] = get_frq_min_max_idx (frq_vec, min_frq_vec, max_frq_vec)
    idx_vec = find((frq_vec >= min_frq_vec) & (frq_vec <= max_frq_vec));
    if ~isempty(idx_vec)
        head_idx = idx_vec(1, 1);
        tail_idx = idx_vec(end, 1);
    else
        head_idx = 1;
        tail_idx = size(frq_vec, 1);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
