% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;
clc;

[cor_z1_mat, cor_z2_mat, raw_tim_cvec, n_sam, sam_frq] = gen_sig();
% save('sample_001.mat', 'cor_z1_mat', 'cor_z2_mat', 'raw_tim_cvec', 'n_sam', 'sam_frq');

% % input
% load sample_001.mat

mean_z1_vec = sum(cor_z1_mat, 2) ./ n_sam;
mean_z2_vec = sum(cor_z2_mat, 2) ./ n_sam;

mean_z1_z2_mat = mean_z1_vec * mean_z2_vec'
z1_z2_mat      = (cor_z1_mat * cor_z2_mat') ./ n_sam
cov_z1_z2_mat  = ((cor_z1_mat - mean_z1_vec) * (cor_z2_mat - mean_z2_vec)') ./ n_sam

cor_r1_mat = real(cor_z1_mat);
cor_m1_mat = imag(cor_z1_mat);
cor_r2_mat = real(cor_z2_mat);
cor_m2_mat = imag(cor_z2_mat);

mean_r1_vec = sum(cor_r1_mat, 2) ./ n_sam;
mean_m1_vec = sum(cor_m1_mat, 2) ./ n_sam;
mean_r2_vec = sum(cor_r2_mat, 2) ./ n_sam;
mean_m2_vec = sum(cor_m2_mat, 2) ./ n_sam;

cov_r1_r2_mat    = ((cor_r1_mat - mean_r1_vec) * (cor_r2_mat - mean_r2_vec).') ./ n_sam;
cov_m1_m2_mat    = ((cor_m1_mat - mean_m1_vec) * (cor_m2_mat - mean_m2_vec).') ./ n_sam;
cov_r1_m2_mat    = ((cor_r1_mat - mean_r1_vec) * (cor_m2_mat - mean_m2_vec).') ./ n_sam;
cov_m1_r2_mat    = ((cor_m1_mat - mean_m1_vec) * (cor_r2_mat - mean_r2_vec).') ./ n_sam;
re_cov_z1_z2_mat = cov_r1_r2_mat + cov_m1_m2_mat
im_cov_z1_z2_mat = cov_r1_m2_mat - cov_m1_r2_mat

% time-domain plot
figure;
subplot(2, 2, 1);
plot_re_time(cor_z1_mat(1, :), raw_tim_cvec);
subplot(2, 2, 2);
plot_re_time(cor_z1_mat(2, :), raw_tim_cvec);
subplot(2, 2, 3);
plot_re_time(cor_z2_mat(1, :), raw_tim_cvec);
subplot(2, 2, 4);
plot_re_time(cor_z2_mat(2, :), raw_tim_cvec);
figure;
subplot(2, 2, 1);
plot_im_time(cor_z1_mat(1, :), raw_tim_cvec);
subplot(2, 2, 2);
plot_im_time(cor_z1_mat(2, :), raw_tim_cvec);
subplot(2, 2, 3);
plot_im_time(cor_z2_mat(1, :), raw_tim_cvec);
subplot(2, 2, 4);
plot_im_time(cor_z2_mat(2, :), raw_tim_cvec);

% histogram itself
figure;
subplot(2, 2, 1);
plot_clpx_hist_pcolr(cor_z1_mat(1, :).', 50, 5);
subplot(2, 2, 2);
plot_clpx_hist_pcolr(cor_z1_mat(2, :).', 50, 5);
subplot(2, 2, 3);
plot_clpx_hist_pcolr(cor_z2_mat(1, :).', 50, 5);
subplot(2, 2, 4);
plot_clpx_hist_pcolr(cor_z2_mat(2, :).', 50, 5);

% histogram with others
figure;
subplot(2, 2, 1);
plot_hist_pcolr(real(cor_z1_mat(1, :).'), real(cor_z1_mat(1, :).'), 50, 5);
subplot(2, 2, 2);
plot_hist_pcolr(real(cor_z1_mat(1, :).'), real(cor_z1_mat(2, :).'), 50, 5);
subplot(2, 2, 3);
plot_hist_pcolr(real(cor_z1_mat(1, :).'), real(cor_z2_mat(1, :).'), 50, 5);
subplot(2, 2, 4);
plot_hist_pcolr(real(cor_z1_mat(1, :).'), real(cor_z2_mat(2, :).'), 50, 5);

% histogram with others
figure;
subplot(2, 2, 1);
plot_hist_pcolr(imag(cor_z1_mat(1, :).'), imag(cor_z1_mat(1, :).'), 50, 5);
subplot(2, 2, 2);
plot_hist_pcolr(imag(cor_z1_mat(1, :).'), imag(cor_z1_mat(2, :).'), 50, 5);
subplot(2, 2, 3);
plot_hist_pcolr(imag(cor_z1_mat(1, :).'), imag(cor_z2_mat(1, :).'), 50, 5);
subplot(2, 2, 4);
plot_hist_pcolr(imag(cor_z1_mat(1, :).'), imag(cor_z2_mat(2, :).'), 50, 5);

% fast fourier transform
figure;
subplot(2, 2, 1);
plot_dft(cor_z1_mat(1, :), sam_frq, 1e6);
subplot(2, 2, 2);
plot_dft(cor_z1_mat(2, :), sam_frq, 1e6);
subplot(2, 2, 3);
plot_dft(cor_z2_mat(1, :), sam_frq, 1e6);
subplot(2, 2, 4);
plot_dft(cor_z2_mat(2, :), sam_frq, 1e6);

function plot_re_time (z_vec, tim_vec)
    plot(tim_vec, real(z_vec));
    
    title('Time-domain Signal');
    xlabel('Time (s)');
    ylabel('Real Number');
    grid on;
end

function plot_im_time (z_vec, tim_vec)
    plot(tim_vec, imag(z_vec));
    
    title('Time-domain Signal');
    xlabel('Time (s)');
    ylabel('Imaginary Number');
    grid on;
end

function plot_hist_pcolr (x_vec, y_vec, n_bin, lim_val)
    
    n_hist_mat = hist3([x_vec, y_vec], 'Nbins', [n_bin, n_bin]);
    n_hist_mat = n_hist_mat.';
    
    n_hist_mat(size(n_hist_mat,1) + 1, size(n_hist_mat, 2) + 1) = 0;
    x_p_vec = linspace(min(x_vec), max(x_vec), size(n_hist_mat, 2));
    y_p_vec = linspace(min(y_vec), max(y_vec), size(n_hist_mat, 1));
    
    pcolor(x_p_vec, y_p_vec, n_hist_mat);
    colormap('hot');
    colorbar;
    xlim([-lim_val, +lim_val]);
    ylim([-lim_val, +lim_val]);
    
    title('Histogram of Complex Variable');
    xlabel('Real Number');
    ylabel('Imaginary Number');
    grid on;
end

function plot_clpx_hist_pcolr (z_vec, n_bin, lim_val)
    x_vec = real(z_vec);
    y_vec = imag(z_vec);
    
    n_hist_mat = hist3([x_vec, y_vec], 'Nbins', [n_bin, n_bin]);
    n_hist_mat = n_hist_mat.';
    
    n_hist_mat(size(n_hist_mat,1) + 1, size(n_hist_mat, 2) + 1) = 0;
    x_p_vec = linspace(min(x_vec), max(x_vec), size(n_hist_mat, 2));
    y_p_vec = linspace(min(y_vec), max(y_vec), size(n_hist_mat, 1));
    
    pcolor(x_p_vec, y_p_vec, n_hist_mat);
    colormap('hot');
    colorbar;
    xlim([-lim_val, +lim_val]);
    ylim([-lim_val, +lim_val]);
    
    title('Histogram of Complex Variable');
    xlabel('Real Number');
    ylabel('Imaginary Number');
    grid on;
end

function plot_dft (z_vec, sam_frq, lim_val)
    fft_vec = fft(z_vec);
    plot((0: length(fft_vec) - 1) * sam_frq / length(fft_vec), abs(fft_vec));
    xlim([0, lim_val]);
    set(gca, 'XScale', 'log');
    
    title('Discrete Fourier Transform (DFT)');
    xlabel('Frequency (Hz)');
    ylabel('DFT Value');
end

function [z1_mat, z2_mat, raw_tim_cvec, n_sam, sam_frq] = gen_sig
    % environment setting
    sam_frq      = 1e6;
    sam_tim      = 1 / sam_frq;
    sim_max_tim  = 100e-3;
    raw_tim_cvec = (0: sam_tim: sim_max_tim);
    n_sam        = size(raw_tim_cvec, 2);
    
    % correlated random coefficient model
%     cor_coef_cvec = 1.5 .* (randn(1, n_sam) + randn(1, n_sam) .* 1i);
    cor_coef_cvec = 1.5 .* randn(1, n_sam) .* (-1 + 1i);
    
    % complex random variable with zero mean
    z1_mat             = zeros(n_sam, 0);
    z1_mat(:, end + 1) = cor_coef_cvec + exp(2 * pi * 1i * 1e3 .* raw_tim_cvec) + (1.00 .* (randn(1, n_sam) + (randn(1, n_sam) .* 1i)));
    z1_mat(:, end + 1) = 2.00 .* (randn(1, n_sam) + (randn(1, n_sam) .* 1i));
    z2_mat             = zeros(n_sam, 0);
    z2_mat(:, end + 1) = cor_coef_cvec + exp(2 * pi * 1i * 1e4 .* raw_tim_cvec) + (1.00 .* (randn(1, n_sam) + (randn(1, n_sam) .* 1i)));
    z2_mat(:, end + 1) = 2.00 .* (randn(1, n_sam) + (randn(1, n_sam) .* 1i));
    
    % number of sample
    n_sam_z1 = size(z1_mat, 2);
    n_sam_z2 = size(z2_mat, 2);
    
    % zeros means assumption
    z1_mat = z1_mat - mean(z1_mat, 1) + 1e-8 .* (rand(n_sam, n_sam_z1) + rand(n_sam, n_sam_z1) .* 1i);
    z2_mat = z2_mat - mean(z2_mat, 1) + 1e-8 .* (rand(n_sam, n_sam_z2) + rand(n_sam, n_sam_z2) .* 1i);
    
    % output
    z1_mat = z1_mat.';
    z2_mat = z2_mat.';
end
