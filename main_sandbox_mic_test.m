% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('capture/');
addpath('miscellaneous/');
addpath('model/');
vic_path = 'output/record/d9-ch8-src4';

% mic index
% mic_src_idx_cvec = [1, 2];
% mic_src_idx_cvec = [1, 2, 3, 6];
% mic_src_idx_cvec = [1, 2, 3, 4, 6, 7];
mic_src_idx_cvec = [1, 2, 3, 4, 5, 6, 7, 8];

% doa configuration
n_src = 4;

% overall configuration
n_mic     = size(mic_src_idx_cvec, 2) / 2;
c         = 343;
d_mic     = 0.09;
lamb      = 2 * d_mic;
cen_frq   = c / lamb;


% other configuration
ang_rad_vec      = (0: 1: 180).' .* pi / 180;
stft_min_frq     = 200;
stft_max_frq     = cen_frq;
ss_stft_min_frq  = 100;
ss_stft_max_frq  = 4000;

% mic position
[~, pos_mat] = l_fwd_ster_fcn([n_mic; 0; n_mic - 1], d_mic, [ndef_ang, ndef_ang], lamb, false);

% wave capture & device configuration
cap_tim          = 10;
dev_name_str     = 'OCTA-CAPTURE';
divr_conf_str    = 'ASIO'; % ASIO4ALL ?
n_chan_val       = 8;
sam_frq_val      = 96e3;
sam_per_fram_val = 256;
bit_dept_str     = '24-bit integer';
mic_map_idx_cvec = [5, 1, 2, 3, 4, 6, 7, 8];

% stft configuration
sam_frq      = sam_frq_val;
win_fcn      = @blackman;
win_size     = fix(50e-3 / (1 / sam_frq));
n_fft        = (2 ^ 14);
hop_size     = ceil(win_size / (2 ^ 4));

% capture
wav_cap_obj = wav_cap_def( ...
    dev_name_str, ...
    divr_conf_str, ...
    n_chan_val, ...
    sam_frq_val, ...
    sam_per_fram_val, ...
    bit_dept_str, ...
    mic_map_idx_cvec, ...
    cap_tim);
wav_cap_obj.init();
wav_cap_obj.cap();
wav_cap_obj.de_init();
src_mat = wav_cap_obj.wav_mat;
wav_cap_obj.save(vic_path);
disp('capture voice');

% select mic
src_mat = src_mat(:, mic_src_idx_cvec);

% sfft procedure
tic;
[stft_mat3, frq_vec, tim_vec]   = stft_mat(src_mat, sam_frq, win_fcn, win_size, n_fft, hop_size);
[post_stft_mat3, post_frq_vec]  = stft_hard_fitr(stft_mat3, frq_vec, stft_min_frq, stft_max_frq);
post_stft_mat3                  = post_stft_mat3 ...
                               ./ (win_size .* (sum(feval(win_fcn, win_size)) / win_size));
toc;
disp('stft');

% IMUSIC's doa procedure (viral!)
cd('doa/1983 - 12 - G. Su/');
tic;
[calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, ang_rad_vec);
toc;
ang_x_deg_vec            = calc_x_rad_vec * 180 / pi;
ang_z_deg_vec            = calc_z_rad_vec * 180 / pi;
[ang_x_deg_vec, idx_vec] = sort(ang_x_deg_vec, 'ascend');
ang_z_deg_vec            = ang_z_deg_vec(idx_vec, 1);
el_deg_vec               = to_el_via_ang_z_deg(ang_z_deg_vec);
az_deg_vec               = to_az_via_ang_x_deg(ang_x_deg_vec, el_deg_vec);
[ang_x_deg_vec, ang_z_deg_vec]
[az_deg_vec, el_deg_vec]
disp('1D IMUSIC');
% tic;
% [calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
% toc;
cd('../..');
% ang_x_deg_vec            = calc_x_rad_vec * 180 / pi;
% ang_z_deg_vec            = calc_z_rad_vec * 180 / pi;
% [ang_x_deg_vec, idx_vec] = sort(ang_x_deg_vec, 'ascend');
% ang_z_deg_vec            = ang_z_deg_vec(idx_vec, 1);
% el_deg_vec               = to_el_via_ang_z_deg(ang_z_deg_vec);
% az_deg_vec               = to_az_via_ang_x_deg(ang_x_deg_vec, el_deg_vec);
% disp('2D IMUSIC');

% doa procedure
cd('doa/now - B. Suksiri/');
tic;
[calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "esprit-auto-pair-x");
toc;
cd('../..');
ang_x_deg_vec            = calc_x_rad_vec * 180 / pi;
ang_z_deg_vec            = calc_z_rad_vec * 180 / pi;
[ang_x_deg_vec, idx_vec] = sort(ang_x_deg_vec, 'ascend');
ang_z_deg_vec            = ang_z_deg_vec(idx_vec, 1);
el_deg_vec               = to_el_via_ang_z_deg(ang_z_deg_vec);
az_deg_vec               = to_az_via_ang_x_deg(ang_x_deg_vec, el_deg_vec);
[ang_x_deg_vec, ang_z_deg_vec]
[az_deg_vec, el_deg_vec]
disp('Proposed Method with ESPRIT');

% doa procedure
cd('doa/now - B. Suksiri/');
tic;
[calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
toc;
cd('../..');
ang_x_deg_vec            = calc_x_rad_vec * 180 / pi;
ang_z_deg_vec            = calc_z_rad_vec * 180 / pi;
[ang_x_deg_vec, idx_vec] = sort(ang_x_deg_vec, 'ascend');
ang_z_deg_vec            = ang_z_deg_vec(idx_vec, 1);
el_deg_vec               = to_el_via_ang_z_deg(ang_z_deg_vec);
az_deg_vec               = to_az_via_ang_x_deg(ang_x_deg_vec, el_deg_vec);
[ang_x_deg_vec, ang_z_deg_vec]
[az_deg_vec, el_deg_vec]
disp('Proposed Method with 2D MUSIC');

% ss test procedure
if false
    % angles
    ang_rad_mat = [calc_x_rad_vec, ndef_ang .* ones(n_src, 1), calc_z_rad_vec];
    
    % frequency selection
    [ss_stft_head_idx, ss_stft_tail_idx] = get_frq_min_max_idx(frq_vec, ss_stft_min_frq, ss_stft_max_frq);
    ss_post_frq_vec                      = frq_vec(ss_stft_head_idx: ss_stft_tail_idx, 1);
    
    % obtain source matrix
    [n_frq_fft, n_tim_fft, n_sam] = size(stft_mat3);
    ss_sig_mat3                   = zeros(n_frq_fft, n_tim_fft, n_src);
    for i_frq = ss_stft_head_idx: 1: ss_stft_tail_idx
        tmp_src_frq   = ss_post_frq_vec(i_frq - ss_stft_head_idx + 1, 1);
        calc_a_fi_mat = l_fwd_ster_frq_fcn( ...
            [n_mic; 0; n_mic - 1], d_mic, ang_rad_mat, tmp_src_frq, cen_frq, lamb, true);
        ss_sig_mat3(i_frq, :, :) = lsqminnorm(calc_a_fi_mat, squeeze(stft_mat3(i_frq, :, :)).').';
    end
    
    % variances of sources
    var_ss_sig_mat = zeros(n_src, n_frq_fft);
    for i_src = 1: 1: n_src
        for i_frq_fft = 1: 1: n_frq_fft
            var_ss_sig_mat(i_src, i_frq_fft) = var(ss_sig_mat3(i_frq_fft, :, i_src));
        end
    end
    figure;
    for i_src = 1: 1: n_src
        plot(frq_vec, var_ss_sig_mat(i_src, :), 'LineWidth', 2);
        hold on;
    end
    
    % isfft procedure
    [syn_mat, ss_tim_vec] = istft_mat(ss_sig_mat3, sam_frq, win_fcn, n_fft, hop_size);
    
    % plot
    figure;
    for i_src = 1: 1: n_src
        plot(ss_tim_vec, syn_mat(:, i_src));
        hold on;
    end
    for i_src = 1: 1: n_src
        soundsc(syn_mat(:, i_src), sam_frq);
        pause(ss_tim_vec(end, 1));
    end
    soundsc(src_mat(:, 1), sam_frq);
end
