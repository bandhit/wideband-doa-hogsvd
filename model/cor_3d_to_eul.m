%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Azimuth and Elevation Conversion
%
% Description : Azimuth and Elevation Conversion
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 18 May 2018, Bandhit Suksiri,
%               Updated: 18 May 2018, Bandhit Suksiri.
%
% Copyright 2018,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [az_deg_val, el_deg_val] = cor_3d_to_eul (cor_a_vec, cor_b_vec)
    x_val      = cor_b_vec(1, 1) - cor_a_vec(1, 1);
    y_val      = cor_b_vec(2, 1) - cor_a_vec(2, 1);
    z_val      = cor_b_vec(3, 1) - cor_a_vec(3, 1);
    r_val      = sqrt((x_val ^ 2) + (y_val ^ 2) + (z_val ^ 2));
    p_val      = sqrt((r_val ^ 2) - (z_val ^ 2));
    az_deg_val = atan2d(y_val, x_val);
    el_deg_val = atan2d(p_val, z_val);

    if az_deg_val < 0
        az_deg_val = 360 + az_deg_val;
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
