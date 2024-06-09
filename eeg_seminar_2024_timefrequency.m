clear all;

% PATH VARS
PATH_EEGLAB = '/home/plkn/eeglab2023.1/';
PATH_AUTOCLEANED = '/mnt/data_fast/sysvself/2_autocleaned/';

% A nice subject list
subject_list = {'Vp02', 'Vp03', 'Vp04', 'Vp05', 'Vp06', 'VP07'};

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% Set complex Morlet wavelet parameters
n_frq = 30;
frqrange = [2, 25];
tfres_range = [400, 100];
EEG = pop_loadset('filename', ['vp_', num2str(2), '_cleaned_tf.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

% Set wavelet time
wtime = -2 : 1 / EEG.srate : 2;

% Determine fft frqs
hz = linspace(0, EEG.srate, length(wtime));

% Create wavelet frequencies and tapering Gaussian widths in temporal domain
tf_freqs = linspace(frqrange(1), frqrange(2), n_frq);
fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

% Init matrices for wavelets
cmw = zeros(length(tf_freqs), length(wtime));
cmwX = zeros(length(tf_freqs), length(wtime));
tlim = zeros(1, length(tf_freqs));

% These will contain the wavelet widths as full width at 
% half maximum in the temporal and spectral domain
obs_fwhmT = zeros(1, length(tf_freqs));
obs_fwhmF = zeros(1, length(tf_freqs));

% Create the wavelets
for frq = 1 : length(tf_freqs)

    % Create wavelet with tapering gaussian corresponding to desired width in temporal domain
    cmw(frq, :) = exp(2 * 1i * pi * tf_freqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);

    % Normalize wavelet
    cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));

    % Create normalized freq domain wavelet
    cmwX(frq, :) = fft(cmw(frq, :)) ./ max(fft(cmw(frq, :)));

    % Determine observed fwhmT
    midt = dsearchn(wtime', 0);
    cmw_amp = abs(cmw(frq, :)) ./ max(abs(cmw(frq, :))); % Normalize cmw amplitude
    obs_fwhmT(frq) = wtime(midt - 1 + dsearchn(cmw_amp(midt : end)', 0.5)) - wtime(dsearchn(cmw_amp(1 : midt)', 0.5));

    % Determine observed fwhmF
    idx = dsearchn(hz', tf_freqs(frq));
    cmwx_amp = abs(cmwX(frq, :)); 
    obs_fwhmF(frq) = hz(idx - 1 + dsearchn(cmwx_amp(idx : end)', 0.5) - dsearchn(cmwx_amp(1 : idx)', 0.5));

end

% Define time window of analysis
pruned_segs = [-500, 2300];
EEG = pop_loadset('filename', ['vp_', num2str(2), '_cleaned_tf.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');
tf_times = EEG.times(dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)));

% ersp matrix
ersps = [];

% Loop subjects
for s = 1 : length(subject_list)

    % Get id
    subject_id = str2num(subject_list{s}(3 : 4));

    % Load data
    EEG = pop_loadset('filename', ['vp_', num2str(subject_id), '_cleaned_tf.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');

    % Get condition indices
    idx_bl  = EEG.trialinfo(:, 4) == 1 & EEG.trialinfo(:, 3) > 0;
    idx_sys = EEG.trialinfo(:, 4) == 2 & EEG.trialinfo(:, 3) > 0;
    idx_slf = EEG.trialinfo(:, 4) == 3 & EEG.trialinfo(:, 3) > 0;

    d = double(EEG.data);

    % tf decomp
    for ch = 1 : size(d, 1)

        % Talk
        fprintf('\ntf decomp subject %i/%i | chan %i/%i...\n', s, numel(subject_list), ch, size(d, 1));

        % Pick channel data
        dch = squeeze(d(ch, :, :));

        % Set convolution length
        convlen = size(dch, 1) * size(dch, 2) + size(cmw, 2) - 1;

        % cmw to freq domain and scale
        cmwX = zeros(n_frq, convlen);
        for f = 1 : n_frq
            cmwX(f, :) = fft(cmw(f, :), convlen);
            cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
        end

        % Get TF-power
        powcube = NaN(n_frq, size(dch, 1), size(dch, 2));
        tmp = fft(reshape(double(dch), 1, []), convlen);
        for f = 1 : n_frq
            as = ifft(cmwX(f, :) .* tmp); 
            as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
            as = reshape(as, size(dch, 1), size(dch, 2));
            powcube(f, :, :) = abs(as) .^ 2;          
        end

        % Cut edge artifacts
        powcube = powcube(:, dsearchn(EEG.times', pruned_segs(1)) : dsearchn(EEG.times', pruned_segs(2)), :);

        % Get condition general baseline values
        ersp_bl = [-500, -200];
        tmp = squeeze(mean(powcube, 3));
        [~, blidx1] = min(abs(tf_times - ersp_bl(1)));
        [~, blidx2] = min(abs(tf_times - ersp_bl(2)));
        blvals = squeeze(mean(tmp(:, blidx1 : blidx2), 2));

        % Calculate ersp
        ersps(s, 1, ch, :, :) = 10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, idx_bl), 3)), blvals));
        ersps(s, 2, ch, :, :) = 10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, idx_sys), 3)), blvals));
        ersps(s, 3, ch, :, :) = 10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, idx_slf), 3)), blvals));

    end % End chanit

end % End subject loop

figure()

subplot(4, 3, 1)
pd = squeeze(mean(squeeze(ersps(:, 1, 15, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Fz - baseline')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 2)
pd = squeeze(mean(squeeze(ersps(:, 2, 15, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Fz - sys')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 3)
pd = squeeze(mean(squeeze(ersps(:, 3, 15, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Fz - self')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 4)
pd = squeeze(mean(squeeze(ersps(:, 1, 65, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('FCz - baseline')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 5)
pd = squeeze(mean(squeeze(ersps(:, 2, 65, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('FCz - sys')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 6)
pd = squeeze(mean(squeeze(ersps(:, 3, 65, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('FCz - self')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 7)
pd = squeeze(mean(squeeze(ersps(:, 1, 16, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Cz - baseline')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 8)
pd = squeeze(mean(squeeze(ersps(:, 2, 16, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Cz - sys')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 9)
pd = squeeze(mean(squeeze(ersps(:, 3, 16, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('Cz - self')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 10)
pd = squeeze(mean(squeeze(ersps(:, 1, 64, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('POz - baseline')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 11)
pd = squeeze(mean(squeeze(ersps(:, 2, 64, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('POz - sys')
clim([-4, 4])
colormap('jet')
xline([0, 1800])

subplot(4, 3, 12)
pd = squeeze(mean(squeeze(ersps(:, 3, 64, :, :)), 1));
contourf(tf_times, tf_freqs, pd, 40, 'linecolor','none')
title('POz - self')
clim([-4, 4])
colormap('jet')
xline([0, 1800])


% Get frequency band indices
idx_theta = tf_freqs >= 4 & tf_freqs <= 7;
idx_alpha = tf_freqs >= 8 & tf_freqs <= 12;
idx_beta = tf_freqs >= 16 & tf_freqs <= 40;

figure()

subplot(3, 3, 1)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 15, idx_theta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 15, idx_theta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 15, idx_theta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('Fz - theta')
xline([0, 1800])

subplot(3, 3, 2)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 15, idx_alpha, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 15, idx_alpha, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 15, idx_alpha, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('Fz - alpha')
xline([0, 1800])

subplot(3, 3, 3)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 15, idx_beta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 15, idx_beta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 15, idx_beta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('Fz - beta')
xline([0, 1800])

subplot(3, 3, 4)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 65, idx_theta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 65, idx_theta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 65, idx_theta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('FCz - theta')
xline([0, 1800])

subplot(3, 3, 5)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 65, idx_alpha, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 65, idx_alpha, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 65, idx_alpha, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('FCz - alpha')
xline([0, 1800])

subplot(3, 3, 6)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 65, idx_beta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 65, idx_beta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 65, idx_beta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('FCz - beta')
xline([0, 1800])

subplot(3, 3, 7)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 64, idx_theta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 64, idx_theta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 64, idx_theta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('POz - theta')
xline([0, 1800])

subplot(3, 3, 8)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 64, idx_alpha, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 64, idx_alpha, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 64, idx_alpha, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('POz - alpha')
xline([0, 1800])

subplot(3, 3, 9)
pd1 = squeeze(mean(squeeze(ersps(:, 1, 64, idx_beta, :)), [1, 2]));
pd2 = squeeze(mean(squeeze(ersps(:, 2, 64, idx_beta, :)), [1, 2]));
pd3 = squeeze(mean(squeeze(ersps(:, 3, 64, idx_beta, :)), [1, 2]));
plot(tf_times, pd1, 'Linewidth', 1.5)
hold on
plot(tf_times, pd2, 'Linewidth', 1.5)
plot(tf_times, pd3, 'Linewidth', 1.5)
legend({'bl', 'sys', 'self'})
title('POz - beta')
xline([0, 1800])
