%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Forward Equation of Steering Vector for L-Shaped Array
%
% Description : Forward Equation of Steering Vector for L-Shaped Array
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 4 June 2017, Bandhit Suksiri,
%               Updated: 4 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ster_mat, pos_mat] = l_fwd_ster_fcn (n_sen_vec, d_sen, ang_rad_mat, lamb, is_eulr_ang)
    norm_frq = 1;
    [ster_mat, pos_mat] = abs_l_fwd_ster_fcn(n_sen_vec, d_sen, ang_rad_mat, ...
                                             norm_frq, norm_frq * lamb, is_eulr_ang);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%