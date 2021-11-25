% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('miscellaneous/');
addpath('simulation/');
addpath('model/');
addpath('doa/1983 - 12 - G. Su/');
addpath('third party/Fast simulation of acoustic room impulse responses/');
revb_path = 'reverberation';
out_path  = 'output';

% load configuration
cfg();

% wave simulation
sim_obj = wav_sim_3d_def();
for info_vec = [az_deg_vec, el_deg_vec, amp_src_vec].'
    az_deg  = info_vec(1, 1);
    el_deg  = info_vec(2, 1);
    amp_src = info_vec(3, 1);
    sim_obj.add_src_ang(az_deg - err_val, el_deg - err_val, amp_src);
end
[~, pos_mat] = l_fwd_ster_fcn(n_mic_vec, d_mic, [ndef_ang, ndef_ang], lamb, false);
for pos_vec = pos_mat.'
    sim_obj.add_pnt(pos_vec);
end
sim_obj.prep();
for tmp_wav_cel = wav_cel
    sim_obj.add_src_wav(tmp_wav_cel{1, 1})
end
if ~isnan(revb_file_name_string)
    load(fullfile(revb_path, revb_file_name_string), 'rir_cel');
    sim_obj.import_rir(rir_cel);
end
sim_obj.sim(sam_frq, sim_max_tim, c);

% measured DOA procedure
meas_x_deg_vec                 = to_ang_x(sim_obj.az_deg_vec * pi / 180, ...
                                          sim_obj.el_deg_vec * pi / 180) * 180 / pi;
meas_z_deg_vec                 = to_ang_z(sim_obj.el_deg_vec * pi / 180) * 180 / pi;
[meas_z_deg_vec, pair_idx_vec] = sort(meas_z_deg_vec, 'ascend');
meas_x_deg_vec                 = meas_x_deg_vec(pair_idx_vec, 1);
n_src                          = size(sim_obj.wav_cel, 1);

% test procedure
n_snr           = size(snr_vec, 1);
calc_x_deg_mat3 = zeros(n_src, n_snr, n_sam);
calc_z_deg_mat3 = zeros(n_src, n_snr, n_sam);
diff_tim_mat    = zeros(n_snr, n_sam);
for i_snr = 1: 1: n_snr
    for i_sam = 1: 1: n_sam % parfor
        % white gaussian noise
        src_mat = sim_obj.sim_wav_mat;
        snr     = snr_vec(i_snr, 1);
        for i_mic = 1: 1: size(src_mat, 2)
            src_mat(:, i_mic) = awgn(src_mat(:, i_mic), snr, 'measured');
        end

        % sfft procedure
        [stft_mat3, frq_vec, tim_vec]   = stft_mat(src_mat, sam_frq, win_fcn, win_size, n_fft, hop_size);
        [post_stft_mat3, post_frq_vec]  = stft_hard_fitr(stft_mat3, frq_vec, stft_min_frq, stft_max_frq);
        post_stft_mat3                  = post_stft_mat3 ...
                                       ./ (win_size .* (sum(feval(win_fcn, win_size)) / win_size));

        % doa procedure
        if is_pair == false
            tim_init = cputime;
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, ang_rad_vec);
            diff_tim = cputime - tim_init;
        else
            tim_init = cputime;
            [calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
            diff_tim = cputime - tim_init
        end
        
        % calculated DOA
        calc_x_deg_vec = calc_x_rad_vec * 180 / pi;
        calc_z_deg_vec = calc_z_rad_vec * 180 / pi;
        
        % save
        calc_x_deg_mat3(:, i_snr, i_sam) = calc_x_deg_vec;
        calc_z_deg_mat3(:, i_snr, i_sam) = calc_z_deg_vec;
        diff_tim_mat(i_snr, i_sam)       = diff_tim;
        
        % percent
        fprintf('>> %6.2f\n', (((i_snr - 1) * n_sam) + i_sam) / (n_snr * n_sam) * 100);
    end
end

% save
save(fullfile(out_path, 'out_1983_12_g_su.mat'), ...
    'is_pair', ...
    'meas_x_deg_vec',  'meas_z_deg_vec', ...
    'calc_x_deg_mat3', 'calc_z_deg_mat3', ...
    'snr_vec', ...
    'diff_tim_mat');
