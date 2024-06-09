clear all;

% PATH VARS
PATH_EEGLAB = '/home/plkn/eeglab2023.1/';
PATH_AUTOCLEANED = '/mnt/data_fast/sysvself/2_autocleaned/';

% A nice subject list
subject_list = {'Vp02', 'Vp03', 'Vp04', 'Vp05', 'Vp06', 'VP07'};

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

% Calculate grand average
erp_ga = squeeze(mean(erps, 1));

% Select channels
channels = [15, 33, 34,...
            65, 19, 20,...
            16,  5,  6,...
            17, 37, 38,...
            ];

figure()
for ch = 1 : length(channels)


    idx_channel = channels(ch);

    pd1 =  squeeze(erp_ga(1, idx_channel, :));
    pd2 =  squeeze(erp_ga(2, idx_channel, :));
    pd3 =  squeeze(erp_ga(3, idx_channel, :));

    subplot(4, 3, ch)
    plot(EEG.times, pd1, 'Linewidth', 1.5)
    hold on
    plot(EEG.times, pd2, 'Linewidth', 1.5)
    plot(EEG.times, pd3, 'Linewidth', 1.5)
    legend({'bl', 'sys', 'self'})
    title(EEG.chanlocs(channels(ch)).labels)
    xline([0, 1800])

end