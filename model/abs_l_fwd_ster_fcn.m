%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Abstract Forward Equation of Steering Vector for L-Shaped Array
%
% Description : Abstract Forward Equation of Steering Vector for L-Shaped Array
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

function [ster_mat, pos_mat] = abs_l_fwd_ster_fcn (n_sen_vec, d_sen, ang_rad_mat, src_frq, c, ...
                                                   is_eulr_ang)
    n_sen_x = n_sen_vec(1, 1);
    n_sen_y = n_sen_vec(2, 1);
    n_sen_z = n_sen_vec(3, 1);
	pos_mat = zeros(1 + n_sen_x + n_sen_y + n_sen_z, 3);
    for i = 1: 1: n_sen_x
        pos_mat(i + 1, :)                     = [d_sen * i, 0, 0];
    end
    for j = 1: 1: n_sen_y
        pos_mat(j + n_sen_x + 1, :)           = [0, d_sen * j, 0];
    end
    for k = 1: 1: n_sen_z
        pos_mat(k + n_sen_y + n_sen_x + 1, :) = [0, 0, d_sen * k];
    end
    ster_mat = abs_fwd_ster_fcn(pos_mat, ang_rad_mat, src_frq, c, is_eulr_ang);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%