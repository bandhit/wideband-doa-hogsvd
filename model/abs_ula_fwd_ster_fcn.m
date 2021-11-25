%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Abstract Forward Equation of Steering Vector for Uniform Linear Array
%
% Description : Abstract Forward Equation of Steering Vector for Uniform Linear Array
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

function [ster_mat, pos_mat] = abs_ula_fwd_ster_fcn (n_sen, d_sen, ang_rad_mat, src_frq, c, ...
                                                     axis_sel)
    is_eulr_ang = true;
    pos_mat     = zeros(n_sen, 3);
    if axis_sel == 'x'
        for i = 1: 1: n_sen
           pos_mat(i, :) = [d_sen * (i - 1), 0, 0];
        end
    elseif axis_sel == 'y'
        for i = 1: 1: n_sen
           pos_mat(i, :) = [0, d_sen * (i - 1), 0];
        end
    elseif axis_sel == 'z'
        for i = 1: 1: n_sen
           pos_mat(i, :) = [0, 0, d_sen * (i - 1)];
        end
    else
        error('Input must be x, y or z.');
    end
    ster_mat = abs_fwd_ster_fcn(pos_mat, ang_rad_mat, src_frq, c, is_eulr_ang);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%