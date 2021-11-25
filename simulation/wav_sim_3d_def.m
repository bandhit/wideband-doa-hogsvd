%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Class of 3-dimensional Geometry Traveling Wave Simulation
%
% Description : 3-dimensional Geometry Traveling Wave Simulation
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 18 May 2017, Bandhit Suksiri,
%               Updated: 24 May 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef wav_sim_3d_def < handle
    properties (Constant = true)
        DIM_SIZE = 3;
    end
	properties (GetAccess = public, SetAccess = protected)
        az_deg_vec    = [];
        el_deg_vec    = [];
        d_src_vec     = [];
        cen_pos_mat   = [];
        norm_mat      = [];
        pnt_mat       = [];
        orth_pnt_mat3 = [];
        wav_cel       = {};
        sim_wav_mat   = [];
        rir_cel       = {};
        sim_wav_mat3  = [];
    end
    methods (Access = public)
        function this = wav_sim_3d_def ()
        end
        function add_src_ang (this, az_deg, el_deg, d_src)
            this.az_deg_vec = [this.az_deg_vec; az_deg];
            this.el_deg_vec = [this.el_deg_vec; el_deg];
            this.d_src_vec  = [this.d_src_vec;  d_src];
        end
        function add_src_wav (this, wav_vec)
            this.wav_cel{end + 1, 1} = wav_vec;
        end
        function add_pnt (this, pnt_vec)
            this.pnt_mat = [this.pnt_mat, pnt_vec];
        end
        function prep (this)
            opts  = optimoptions(@fsolve, 'Display','off', 'Algorithm', 'levenberg-marquardt');
            n_src = size(this.d_src_vec, 1);
            n_pnt = size(this.pnt_mat,   2);
            this.cen_pos_mat   = zeros(this.DIM_SIZE, n_src);
            this.orth_pnt_mat3 = zeros(this.DIM_SIZE, n_src, n_pnt);
            this.rir_cel       = cell(n_src, n_pnt);
            for i = 1: 1: n_src
                az_rad = this.az_deg_vec(i, 1) * pi / 180;
                el_rad = this.el_deg_vec(i, 1) * pi / 180;
                d_src  = this.d_src_vec(i, 1);
                cen_pos_vec = d_src .* [sin(el_rad) * cos(az_rad);
                                        sin(el_rad) * sin(az_rad);
                                        cos(el_rad)];
                cen_pos_vec(abs(cen_pos_vec) <= eps) = 0;
                norm_vec = - cen_pos_vec ./ norm(cen_pos_vec, 2);
                norm_vec(abs(norm_vec) <= eps) = 0;
                sys_fcn  = @(orth_pnt_vec) ((norm_vec' * cen_pos_vec) - ...
                                            (norm_vec' * orth_pnt_vec));

                for j = 1: 1: n_pnt
                    pnt_vec         = this.pnt_mat(:, j);
                    orth_fcn        = @(orth_pnt_vec) ((pnt_vec - orth_pnt_vec)' * ...
                                                       (cen_pos_vec - orth_pnt_vec));
                    norm_x_orth_fcn = @(orth_pnt_vec) ((pnt_vec(1, 1) - orth_pnt_vec(1, 1)) - ...
                                                      (norm(pnt_vec - orth_pnt_vec, 2) .* ...
                                                       norm_vec(1, 1)));
                    norm_y_orth_fcn = @(orth_pnt_vec) ((pnt_vec(2, 1) - orth_pnt_vec(2, 1)) - ...
                                                      (norm(pnt_vec - orth_pnt_vec, 2) .* ...
                                                       norm_vec(2, 1)));
                    norm_z_orth_fcn = @(orth_pnt_vec) ((pnt_vec(3, 1) - orth_pnt_vec(3, 1)) - ...
                                                      (norm(pnt_vec - orth_pnt_vec, 2) .* ...
                                                       norm_vec(3, 1)));
                    obj_fcn         = @(orth_pnt_vec) ([sys_fcn(orth_pnt_vec);
                                                        orth_fcn(orth_pnt_vec);
                                                        norm_x_orth_fcn(orth_pnt_vec);
                                                        norm_y_orth_fcn(orth_pnt_vec);
                                                        norm_z_orth_fcn(orth_pnt_vec)]);
                    orth_pnt_vec    = fsolve(obj_fcn, cen_pos_vec, opts);
                    
                    orth_pnt_vec(abs(orth_pnt_vec) <= eps) = 0;
                    this.orth_pnt_mat3(:, i, j) = orth_pnt_vec;
                end
                
                this.cen_pos_mat(:, i) = cen_pos_vec;
                this.norm_mat(:, i)    = norm_vec;
                this.rir_cel{i, j}     = NaN;
            end
            
        end
        function import_rir (this, rir_cel)
            this.rir_cel = rir_cel;
        end
        function sim (this, sam_frq, max_tim, sped)
            sam_tim           = 1 / sam_frq;
            n_src             = size(this.d_src_vec, 1);
            n_pnt             = size(this.pnt_mat,   2);
            n_tim             = fix(max_tim / sam_tim);
            this.sim_wav_mat  = zeros(n_tim, n_pnt);
            this.sim_wav_mat3 = zeros(n_tim, n_pnt, n_src);
            for i = 1: 1: n_src
                each_sim_wav_mat = zeros(n_tim, n_pnt);
                org_wav_vec      = this.wav_cel{i, 1};
                for j = 1: 1: n_pnt
                    rir_obj = this.rir_cel{i, j};
                    if ~isnan(rir_obj)
                        rir_vec           = rir_obj;
                        norm_rir_vec      = rir_vec ./ max(abs(rir_vec));
                        [~, idx_rir_max]  = max(norm_rir_vec);
                        shif_norm_rir_vec = norm_rir_vec(idx_rir_max: end, 1);
                        wav_vec = conv(shif_norm_rir_vec, org_wav_vec);
                    else
                        wav_vec = org_wav_vec;
                    end
                    wav_size = size(wav_vec, 1);
                    
                    pnt_vec        = this.pnt_mat(:, j);
                    orth_pnt_vec   = this.orth_pnt_mat3(:, i, j);
                    d_orth_pnt     = norm(orth_pnt_vec - pnt_vec, 2);
                    tim_orth_pnt   = d_orth_pnt / sped;
                    n_tim_orth_pnt = fix(tim_orth_pnt / sam_tim);
                    if (wav_size + n_tim_orth_pnt) <= n_tim
                        each_sim_wav_mat((1: wav_size) + n_tim_orth_pnt, j) = ...
                        each_sim_wav_mat((1: wav_size) + n_tim_orth_pnt, j) + ...
                            wav_vec;
                    elseif n_tim_orth_pnt < n_tim
                        each_sim_wav_mat(n_tim_orth_pnt + 1: n_tim, j) = ...
                        each_sim_wav_mat(n_tim_orth_pnt + 1: n_tim, j) + ...
                            wav_vec(1: n_tim - n_tim_orth_pnt, 1);
                    end
                end
                this.sim_wav_mat           = this.sim_wav_mat + each_sim_wav_mat;
                this.sim_wav_mat3(:, :, i) = each_sim_wav_mat;
            end
        end
        function plot (this)
            n_src   = size(this.d_src_vec, 1);
            n_pnt   = size(this.pnt_mat,   2);
            clr_mat = hsv(n_src);
            for i = 1: 1: n_src
                clr_col     = clr_mat(i, :);
                cen_pos_vec = this.cen_pos_mat(:, i);
                plot3([0; cen_pos_vec(1, 1)], ...
                      [0; cen_pos_vec(2, 1)], ...
                      [0; cen_pos_vec(3, 1)], ...
                      '--o', ...
                      'Color',           clr_col, ...
                      'MarkerFaceColor', clr_col);
                hold on;
                for j = 1: 1: n_pnt
                    pnt_vec      = this.pnt_mat(:, j);
                    orth_pnt_vec = this.orth_pnt_mat3(:, i, j);
                    plot3([pnt_vec(1, 1); orth_pnt_vec(1, 1)], ...
                          [pnt_vec(2, 1); orth_pnt_vec(2, 1)], ...
                          [pnt_vec(3, 1); orth_pnt_vec(3, 1)], ...
                          '-o', ...
                          'Color',           clr_col, ...
                          'MarkerFaceColor', clr_col);
                    hold on;
                end
                min_x_plne = min(min(this.orth_pnt_mat3(1, i, :)));
                max_x_plne = max(max(this.orth_pnt_mat3(1, i, :)));
                min_y_plne = min(min(this.orth_pnt_mat3(2, i, :)));
                max_y_plne = max(max(this.orth_pnt_mat3(2, i, :)));
                norm_vec   = this.norm_mat(:, i);
                if ~(norm_vec(1, 1) == 0) && ...
                   ~(norm_vec(2, 1) == 0) && ...
                   ~(norm_vec(3, 1) == 0)
                    [x_mat, y_mat] = meshgrid(linspace(min_x_plne - 0.5, max_x_plne + 0.5), ...
                                              linspace(min_y_plne - 0.5, max_y_plne + 0.5), ...
                                              10);
                    z_mat = ((norm_vec' * cen_pos_vec)   - ...
                             (norm_vec(1, 1) .* x_mat)   - ...
                             (norm_vec(2, 1) .* y_mat)) ./ norm_vec(3, 1);
                    clr_mat3          = ones(size(z_mat, 1), size(z_mat, 2), this.DIM_SIZE);
                    clr_mat3(:, :, 1) = clr_col(1, 1) .* clr_mat3(:, :, 1);
                    clr_mat3(:, :, 2) = clr_col(1, 2) .* clr_mat3(:, :, 2);
                    clr_mat3(:, :, 3) = clr_col(1, 3) .* clr_mat3(:, :, 3);
                    surf(x_mat, y_mat, z_mat, clr_mat3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                end
            end
            scatter3(this.pnt_mat(1, :), this.pnt_mat(2, :), this.pnt_mat(3, :), 'black', 'filled');
            grid on;
            grid minor;
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
