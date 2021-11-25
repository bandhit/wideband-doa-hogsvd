%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Eyring's Formula
%
% Description : Eyring's Formula RT-60
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 14 September 2018, Bandhit Suksiri,
%               Updated: 14 September 2018, Bandhit Suksiri.
%
% Copyright 2017 - 2018,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rt_val = eyr_rt_fcn (c, room_dim_vec, wal_ref_vec)
    v  =  room_dim_vec(1, 1) * room_dim_vec(2, 1) * room_dim_vec(3, 1);
    sa = (room_dim_vec(1, 1) * room_dim_vec(3, 1) * wal_ref_vec(1, 1)) + ...
         (room_dim_vec(1, 1) * room_dim_vec(3, 1) * wal_ref_vec(2, 1)) + ...
         (room_dim_vec(2, 1) * room_dim_vec(3, 1) * wal_ref_vec(3, 1)) + ...
         (room_dim_vec(2, 1) * room_dim_vec(3, 1) * wal_ref_vec(4, 1)) + ...
         (room_dim_vec(1, 1) * room_dim_vec(2, 1) * wal_ref_vec(5, 1)) + ...
         (room_dim_vec(1, 1) * room_dim_vec(2, 1) * wal_ref_vec(6, 1));
    s  = (room_dim_vec(1, 1) * room_dim_vec(3, 1) * 2) + ...
         (room_dim_vec(2, 1) * room_dim_vec(3, 1) * 2) + ...
         (room_dim_vec(1, 1) * room_dim_vec(2, 1) * 2);
    rt_val = (24 * log(10) * v) / (c * ( - s * log(1 - (sa / s))));
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%