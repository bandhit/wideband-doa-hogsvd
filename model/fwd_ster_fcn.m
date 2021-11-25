%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Abstract Forward Equation of Steering Vector
%
% Description : Abstract Forward Equation of Steering Vector
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

function ster_mat = fwd_ster_fcn (pos_mat, ang_rad_mat, lamb, is_eulr_ang)
    norm_frq = 1;
    ster_mat = abs_fwd_ster_fcn(pos_mat, ang_rad_mat, norm_frq, norm_frq * lamb, is_eulr_ang);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%