%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Forward Equation of Steering Vector for Uniform Linear Array
%
% Description : Forward Equation of Steering Vector for Uniform Linear Array
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 21 June 2017, Bandhit Suksiri,
%               Updated: 21 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ster_mat, pos_mat] = ula_fwd_ster_frq_fcn (n_sen, d_sen, ang_rad_mat, ...
                                                     src_frq, cen_frq, lamb, axis_sel)
    [ster_mat, pos_mat] = abs_ula_fwd_ster_fcn(n_sen, d_sen, ang_rad_mat, ...
                                               src_frq, cen_frq * lamb, axis_sel);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%