clear all;

% PATH VARS
PATH_EEGLAB = '/home/plkn/eeglab2023.1/';
PATH_RAW = '/mnt/data_fast/sysvself/0_eeg/';
PATH_LOGFILES = '/mnt/data_fast/sysvself/logfiles/';
PATH_ICSET = '/mnt/data_fast/sysvself/1_icset/';
PATH_AUTOCLEANED = '/mnt/data_fast/sysvself/2_autocleaned/';

% A nice subject list
subject_list = {'Vp02', 'Vp03', 'Vp04', 'Vp05', 'Vp06', 'VP07'};
subject_list = {'VP07'};

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% Get chanlocfile
channel_location_file = which('standard-10-5-cap385.elp');

% Loop subjects
for s = 1 : length(subject_list)

    % Load EEG
    EEG = pop_loadbv(PATH_RAW, [subject_list{s}, '.vhdr'], [], []);

    % Get id
    subject_id = str2num(subject_list{s}(3 : 4));

    % Read log
    LOG = readtable([PATH_LOGFILES, 'VP', subject_list{s}(3 : 4), '_logEdgy.txt'], "NumHeaderLines", 3);

    % Iterate events
    trialinfo = [];
    trial_nr_total = 0;
    for e = 1 : length(EEG.event)

        % If baseline 1
        if strcmpi(EEG.event(e).type, 'S 10')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            1];
        end

        % If baseline 2
        if strcmpi(EEG.event(e).type, 'S 20')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            1];
        end

        % If sys 1
        if strcmpi(EEG.event(e).type, 'S 30')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            2];
        end

        % If sys 2
        if strcmpi(EEG.event(e).type, 'S 40')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            2];
        end


        % If self 1
        if strcmpi(EEG.event(e).type, 'S 50')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            3];
        end

        % If self 2
        if strcmpi(EEG.event(e).type, 'S 60')
            EEG.event(e).type = 'X';
            trial_nr_total = trial_nr_total + 1;
            trialinfo(trial_nr_total, :) = [e,...
                                            subject_id,...
                                            trial_nr_total,...
                                            3];
        end
    end

    % Check trialcount
    if trial_nr_total ~= 600
        fprintf('\n\n\nSOMETHING IS WEIIIRDDD with the trials!!!!!!\n\n\n');
        pause;
    end

    % Add to EEG
    EEG.trialinfo = trialinfo;

    % Select EEG channels
    EEG = pop_select(EEG, 'channel', [1 : 64]);

    % Add FCz as empty channel
    EEG.data(end + 1, :) = 0;
    EEG.nbchan = size(EEG.data, 1);
    EEG.chanlocs(end + 1).labels = 'FCz';

    % Add channel locations
    EEG = pop_chanedit(EEG, 'lookup', channel_location_file);

    % Save original channel locations (for later interpolation)
    EEG.chanlocs_original = EEG.chanlocs;

    % Reref to CPz, so that FCz obtains non-interpolated data
    EEG = pop_reref(EEG, 'CPz');

    % Resample data
    EEG_TF = pop_resample(EEG, 200);

    % Filter
    EEG    = pop_basicfilter(EEG,    [1 : EEG.nbchan],    'Cutoff', [0.01, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 6, 'RemoveDC', 'on', 'Boundary', 'boundary'); 
    EEG_TF = pop_basicfilter(EEG_TF, [1 : EEG_TF.nbchan], 'Cutoff', [   1, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 6, 'RemoveDC', 'on', 'Boundary', 'boundary');
        
    % Bad channel detection
    [EEG, EEG.chans_rejected]       = pop_rejchan(EEG,    'elec', [1 : EEG.nbchan],    'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    [EEG_TF, EEG_TF.chans_rejected] = pop_rejchan(EEG_TF, 'elec', [1 : EEG_TF.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');

    % Interpolate channels
    EEG    = pop_interp(EEG,    EEG.chanlocs_original,    'spherical');
    EEG_TF = pop_interp(EEG_TF, EEG_TF.chanlocs_original, 'spherical');

    % Reref common average
    EEG    = pop_reref(EEG,    []);
    EEG_TF = pop_reref(EEG_TF, []);

    % Determine rank of data
    dataRank = sum(eig(cov(double(EEG_TF.data'))) > 1e-6);

    % Epoch EEG data
    [EEG, idx_to_keep] = pop_epoch(EEG, {'X'}, [-0.3, 2.8], 'newname', ['vp_', num2str(subject_id), '_epoched'], 'epochinfo', 'yes');
    EEG.trialinfo =  EEG.trialinfo(idx_to_keep, :);
    EEG = pop_rmbase(EEG, [-200, 0]);
    [EEG_TF, idx_to_keep] = pop_epoch(EEG_TF, {'X'}, [-0.8, 3.3], 'newname', ['vp_', num2str(subject_id), '_epoched'],  'epochinfo', 'yes');
    EEG_TF.trialinfo =  EEG_TF.trialinfo(idx_to_keep, :);
    EEG_TF = pop_rmbase(EEG_TF, [-200, 0]);

    % Autoreject trials in tf-set
    [EEG_TF, EEG_TF.rejected_epochs] = pop_autorej(EEG_TF, 'nogui', 'on');

    % Remove rejected trials from trialinfo of tf-set
    EEG_TF.trialinfo(EEG_TF.rejected_epochs, :) = [];

    % Runica & ICLabel
    EEG_TF = pop_runica(EEG_TF, 'extended', 1, 'interrupt', 'on', 'PCA', dataRank);
    EEG_TF = iclabel(EEG_TF);

    % Find nobrainer
    EEG_TF.nobrainer = find(EEG_TF.etc.ic_classification.ICLabel.classifications(:, 1) < 0.3 | EEG_TF.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3);

    % Copy ICs to erpset
    EEG = pop_editset(EEG, 'icachansind', 'EEG_TF.icachansind', 'icaweights', 'EEG_TF.icaweights', 'icasphere', 'EEG_TF.icasphere');
    EEG.etc = EEG_TF.etc;
    EEG.nobrainer = EEG_TF.nobrainer;

    % Save IC set
    pop_saveset(EEG, 'filename', ['vp_', num2str(subject_id), '_icset_erp.set'], 'filepath', PATH_ICSET, 'check', 'on');
    pop_saveset(EEG_TF, 'filename', ['vp_', num2str(subject_id), '_icset_tf.set'], 'filepath', PATH_ICSET, 'check', 'on');

    % Remove components
    EEG    = pop_subcomp(EEG, EEG.nobrainer, 0);
    EEG_TF = pop_subcomp(EEG_TF, EEG_TF.nobrainer, 0);

    % trial rejection for erp-set
    [EEG, EEG.rejected_epochs] = pop_autorej(EEG, 'nogui', 'on');
    EEG.trialinfo(EEG.rejected_epochs, :) = [];

    % Save clean data
    pop_saveset(EEG, 'filename', ['vp_', num2str(subject_id), '_cleaned_erp.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');
    pop_saveset(EEG_TF, 'filename', ['vp_', num2str(subject_id), '_cleaned_tf.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');

end % End subject loop


