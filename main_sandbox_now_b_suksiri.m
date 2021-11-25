% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;
clc;

% include
addpath('miscellaneous/');
addpath('simulation/');
addpath('model/');
addpath('doa/now - B. Suksiri/');
addpath('source separation/direct/');
addpath('third party/Fast simulation of acoustic room impulse responses/');
addpath('third party/KN Rank Estimation/');
revb_path = 'reverberation';
out_path  = 'output';

% load configuration
cfg;
doa_stft_min_frq = stft_min_frq;
doa_stft_max_frq = stft_max_frq;
% doa_stft_min_frq = 90;
% doa_stft_max_frq = 10010;
ss_stft_min_frq  = 100;
ss_stft_max_frq  = 8000;

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
                                          sim_obj.el_deg_vec * pi / 180) * 180 / pi
meas_z_deg_vec                 = to_ang_z(sim_obj.el_deg_vec * pi / 180) * 180 / pi
[meas_z_deg_vec, pair_idx_vec] = sort(meas_z_deg_vec, 'ascend');
meas_x_deg_vec                 = meas_x_deg_vec(pair_idx_vec, 1);
n_src                          = size(sim_obj.wav_cel, 1);

% doa test procedure
if true
    src_mat = sim_obj.sim_wav_mat;

    snr = 20;
    for i_mic = 1: 1: size(src_mat, 2)
        src_mat(:, i_mic) = awgn(src_mat(:, i_mic), snr, 'measured');
    end

    % sfft procedure
    [stft_mat3, frq_vec, tim_vec] = stft_mat(src_mat, sam_frq, win_fcn, win_size, n_fft, hop_size);
    
    % frequency selection
    [doa_stft_head_idx, doa_stft_tail_idx] = get_frq_min_max_idx(frq_vec, doa_stft_min_frq, doa_stft_max_frq);
    doa_post_frq_vec                       = frq_vec(doa_stft_head_idx: doa_stft_tail_idx, 1);
    doa_post_stft_mat3                     = stft_mat3(doa_stft_head_idx: doa_stft_tail_idx, :, :);
    doa_post_stft_mat3                     = doa_post_stft_mat3 ./ (win_size .* (sum(feval(win_fcn, win_size)) / win_size));
    
    % doa procedure
    tim_init = cputime;
%     [ang_x_rad_vec, ang_z_rad_vec] = doa_frq_n_src(doa_post_stft_mat3, doa_post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "esprit-auto-pair-src-var");
%     [ang_x_rad_vec, ang_z_rad_vec] = doa_frq_n_src(doa_post_stft_mat3, doa_post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "music-auto-pair-norm");
    [ang_x_rad_vec, ang_z_rad_vec] = doa_frq_n_src(doa_post_stft_mat3, doa_post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "music-signal-auto-pair-norm");
    diff_tim = cputime - tim_init
    ang_x_rad_vec * 180 / pi
    ang_z_rad_vec * 180 / pi
end

% ss test procedure
if false
    % angles
    ang_rad_mat = [ang_x_rad_vec, ndef_ang .* ones(n_src, 1), ang_z_rad_vec];
%     ang_rad_mat = [meas_x_deg_vec / 180 * pi, ndef_ang .* ones(n_src, 1), meas_z_deg_vec / 180 * pi];
    
    % frequency selection
    [ss_stft_head_idx, ss_stft_tail_idx] = get_frq_min_max_idx(frq_vec, ss_stft_min_frq, ss_stft_max_frq);
    ss_post_frq_vec                      = frq_vec(ss_stft_head_idx: ss_stft_tail_idx, 1);
    
    % obtain inverse mixture matrix
    inv_mix_mat3 = get_inv_mix_mat(ang_rad_mat, ss_post_frq_vec, n_mic, d_mic, lamb, cen_frq);
    
    % obtain source matrix
    [n_frq_fft, n_tim_fft, n_sam] = size(stft_mat3);
    ss_sig_mat3                   = zeros(n_frq_fft, n_tim_fft, n_src);
    for i_frq = ss_stft_head_idx: 1: ss_stft_tail_idx
        ss_sig_mat3(i_frq, :, :) = (squeeze(inv_mix_mat3(:, :, i_frq - ss_stft_head_idx + 1)) * squeeze(stft_mat3(i_frq, :, :)).').';
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
        pause(sim_max_tim);
        soundsc(syn_mat(:, i_src), sam_frq);
    end
end

% signal from sources
if false
    org_src_mat = squeeze(sim_obj.sim_wav_mat3(:, 1, :));

    % sfft procedure
    org_stft_mat3                           = stft_mat(org_src_mat, sam_frq, win_fcn, win_size, n_fft, hop_size);
    [post_org_stft_mat3, org_post_frq_vec]  = stft_hard_fitr(org_stft_mat3, frq_vec, doa_stft_min_frq, doa_stft_max_frq);
    post_org_stft_mat3                      = post_org_stft_mat3 ...
                                           ./ (win_size .* (sum(feval(win_fcn, win_size)) / win_size));

    % variances of sources
    n_org_post_frq  = size(org_post_frq_vec, 1);
    var_org_src_mat = zeros(n_src, n_org_post_frq);
    for i_src = 1: 1: n_src
        for i_post_frq = 1: 1: n_org_post_frq
            var_org_src_mat(i_src, i_post_frq) = var(post_org_stft_mat3(i_post_frq, :, i_src));
        end
    end
%     figure;
%     for i_src = 1: 1: n_src
%         plot(org_post_frq_vec, var_org_src_mat(i_src, :), 'LineWidth', 2);
%         hold on;
%     end
    
    i_src = 1;
    co_var_mat = zeros(n_org_post_frq, n_org_post_frq);
    for i_frq_1 = 1: 1: n_org_post_frq
        for i_frq_2 = 1: 1: n_org_post_frq
            co_var_mat(i_frq_1, i_frq_2) = (post_org_stft_mat3(i_frq_1, :, i_src) * post_org_stft_mat3(i_frq_2, :, i_src)') ./ n_org_post_frq;
        end
    end
    fig_obj = figure;
    subplot(2, 1, 1);
    plot_co_var (org_post_frq_vec, real(co_var_mat));
    title({'Cross-covariance between a source and itself', 'with distinct frequencies (real part)'});
    subplot(2, 1, 2);
    plot_co_var (org_post_frq_vec, imag(co_var_mat));
    title({'Cross-covariance between a source and itself', 'with distinct frequencies (imaginary part)'});
%     save_png('co_var', fig_obj, 13, 20, '300');
end

function plot_co_var (org_post_frq_vec, co_var_mat)
    abs_co_var_mat = log10(abs(co_var_mat) + 1e-9);
    
    h_obj = pcolor(org_post_frq_vec, org_post_frq_vec, abs_co_var_mat);
    set(gca, 'xscale', 'log', 'yscale', 'log');
    set(h_obj, 'EdgeColor', 'none');
    colormap('hot');
    c_obj = colorbar;
    shading interp;
    caxis([-9, max(max(abs_co_var_mat))]);
    c_obj.Label.String = 'Cross-covariance value (absolute log scale)';
    xlabel('Frequency (Hz)');
    ylabel('Frequency (Hz)');
end

