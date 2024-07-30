clear all;

% PATH VARS
PATH_EEGLAB = '/home/plkn/eeglab2024.0/';
PATH_AUTOCLEANED = '/mnt/data_fast/sysvself/2_autocleaned/';

% A nice subject list
subject_list = {'Vp02', 'Vp03', 'Vp04', 'Vp05', 'Vp06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', 'VP12', 'Vp13'};

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% ERP matrix
erps = [];

% Loop subjects
for s = 1 : length(subject_list)

    % Get id
    subject_id = str2num(subject_list{s}(3 : 4));

    % Load data
    EEG = pop_loadset('filename', ['vp_', num2str(subject_id), '_cleaned_erp.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');

    % Get condition indices
    idx_bl  = EEG.trialinfo(:, 4) == 1 & EEG.trialinfo(:, 3) > 0;
    idx_sys = EEG.trialinfo(:, 4) == 2 & EEG.trialinfo(:, 3) > 0;
    idx_slf = EEG.trialinfo(:, 4) == 3 & EEG.trialinfo(:, 3) > 0;

    % Calculate ERPs
    erps(s, 1, :, :) = squeeze(mean(EEG.data(:, :, idx_bl), 3));
    erps(s, 2, :, :) = squeeze(mean(EEG.data(:, :, idx_sys), 3));
    erps(s, 3, :, :) = squeeze(mean(EEG.data(:, :, idx_slf), 3));

end

% Prune
keep_idx = EEG.times <= 2300;
erps = erps(:, :, :, keep_idx);
erp_times = EEG.times(keep_idx);

% Select channels
chan_idx = [15, 65, 19, 20, 16];

% Average across electrodes
frontal_erps = squeeze(mean(erps(:, :, chan_idx, :), 3));

% Time indices for parametrization
time_idx = erp_times >= 1300 & erp_times < 1800;

% Average across participants
ga_erps =  squeeze(mean(frontal_erps, 1));
ga_topo = squeeze(mean(erps(:, :, :, time_idx), [1, 4]));

% Plot lineplot
figure()
plot(erp_times, ga_erps(1, :), 'k', 'LineWidth', 2.5)
hold on
plot(erp_times, ga_erps(2, :), 'r', 'LineWidth', 2.5)
plot(erp_times, ga_erps(3, :), 'g', 'LineWidth', 2.5)
legend({'baseline', 'system', 'self'})
xlabel('ms')
ylabel('mV')
xline([1300, 1800])
title('frontal ERP (Fz, FCz, FC1, FC2, Cz)')

% Plot topos
figure()
subplot(1, 3, 1)
topoplot(ga_topo(1, :), EEG.chanlocs)
title('baseline')
clim([-4, 4])
subplot(1, 3, 2)
topoplot(ga_topo(2, :), EEG.chanlocs)
title('system')
clim([-4, 4])
subplot(1, 3, 3)
topoplot(ga_topo(3, :), EEG.chanlocs)
title('self')
clim([-4, 4])

% Some stats
res = [];
counter = 0;
for s = 1 : length(subject_list)
    id = subject_list{s};
    id = str2num(id(3 : 4));
    for cond = 1 : 3
        counter = counter + 1;
        amp = mean(squeeze(frontal_erps(s, cond, time_idx)))
        res(counter, :) = [id, cond, amp];
    end
end

% Save
writematrix(res, 'erp_results.csv');











