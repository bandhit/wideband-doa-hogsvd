%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : ESPRIT Angle Function on ULA Model
%
% Description : ESPRIT Angle Function on ULA Model, IEEE
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 14 March 2019, Bandhit Suksiri,
%               Updated: 14 March 2019, Bandhit Suksiri.
%
% Copyright 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_ang_x_rad_vec = esprit_ang_fcn (eig_val_x_vec, fc_src_frq, cen_frq, lamb, d_sen)
    calc_ang_x_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * fc_src_frq)) * angle(eig_val_x_vec));
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%