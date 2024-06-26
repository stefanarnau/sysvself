% This command clears the matlab workspace
clear all;

% Some path variables are defined here.
PATH_EEGLAB = '/home/plkn/eeglab2024.0/';
PATH_RAW = '/mnt/data_fast/sysvself/0_eeg/';
PATH_LOGFILES = '/mnt/data_fast/sysvself/logfiles/';
PATH_ICSET = '/mnt/data_fast/sysvself/1_icset/';
PATH_AUTOCLEANED = '/mnt/data_fast/sysvself/2_autocleaned/';

% A list of file names for eeg data
subject_list = {'Vp02', 'Vp03', 'Vp04', 'Vp05', 'Vp06', 'VP07', 'VP08'};
subject_list = {'VP08'};

% Initialize eeglab and add path for channel location file
addpath(PATH_EEGLAB);
addpath([PATH_EEGLAB, 'plugins/dipfit/standard_BESA/']);
eeglab;

% Get channel location information
channel_locations = which('standard-10-5-cap385.elp');

% Loop datasets
for s = 1 : length(subject_list)

    % Load EEG data from file(s)
    EEG = pop_loadbv(PATH_RAW, [subject_list{s}, '.vhdr'], [], []);

    % an identifier as an integer from string
    subject_id = str2num(subject_list{s}(3 : 4));

    % Read log-file (although we are doing nothing with it so far...)
    LOG = readtable([PATH_LOGFILES, 'VP', subject_list{s}(3 : 4), '_logEdgy.txt'], "NumHeaderLines", 3);

    % Iterate events and build a trialinfo list. This list then contains the experimental condition of
    % each trial. We also rename all events that should be used for time-locking (epoching) to 'X' (could
    % be anything, but Elon would be proud...)
    trialinfo = [];
    trial_nr_total = 0;
    for e = 1 : length(EEG.event)

        % The event numbers 10, 20, 30 40, 50, 60 code for the onset of
        % the memory array (the 5 colored circles). We use those events 
        % for time-locking here. 10 and 20 code for baseline blocks, 30 and 40
        % code for "system" condition blocks (where the response could be flipped), and
        % 50 and 60 code for "self" condition blocks (where the target color could be 
        % flipped).

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

    % Check trialcount. There should be 600 trials. If not, something went wrong 
    if trial_nr_total ~= 600
        error('\nSomething went wrong...\n');
    end

    % We add the trialinfo array to the 'EEG' data structure
    EEG.trialinfo = trialinfo;

    % Add FCz (during recording this was our reference channel) as an empty channel (all zeros)
    EEG.data(end + 1, :) = 0;
    EEG.nbchan = size(EEG.data, 1);
    EEG.chanlocs(end + 1).labels = 'FCz';

    % Add channel locations (coordinates obtained from file)
    EEG = pop_chanedit(EEG, 'lookup', channel_locations);

    % Save original channel locations (for later interpolation, as some channels might be dropped)
    EEG.chanlocs_original = EEG.chanlocs;

    % Rereference the data to CPz, so that FCz obtains non-interpolated data (or any data at all, all zeros so far).
    % This also means that CPz is zero now. 
    EEG = pop_reref(EEG, 'CPz');

    % For time-frequency analysis, we resample the data to 200 Hz. This reduces the amount of data to 1/5 compared to
    % the original 1000 Hz while being still more than sufficient for frequencies up to ~50 Hz.
    EEG_TF = pop_resample(EEG, 200);

    % Filter the data. For ERP analysis we filter the data with a bandpass butterworth-filter from 0.01 Hz to 30 Hz. For time-frequency analysis the filter-boundaries are 1 Hz to 30 Hz
    EEG    = pop_basicfilter(EEG,    [1 : EEG.nbchan],    'Cutoff', [0.01, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 6, 'RemoveDC', 'on', 'Boundary', 'boundary'); 
    EEG_TF = pop_basicfilter(EEG_TF, [1 : EEG_TF.nbchan], 'Cutoff', [   1, 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 6, 'RemoveDC', 'on', 'Boundary', 'boundary');
        
    % Bad channel detection. A sfunction detecting and removing channels with bad data quality based on statistical properties.
    [EEG, EEG.chans_rejected]       = pop_rejchan(EEG,    'elec', [1 : EEG.nbchan],    'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    [EEG_TF, EEG_TF.chans_rejected] = pop_rejchan(EEG_TF, 'elec', [1 : EEG_TF.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');

    % Interpolate missing channels (That is the previously rejected channels)
    EEG    = pop_interp(EEG,    EEG.chanlocs_original,    'spherical');
    EEG_TF = pop_interp(EEG_TF, EEG_TF.chanlocs_original, 'spherical');

    % Rereference the data to common average reference. This means, at each time point, subtracting the average of all channels from each channel.
    EEG    = pop_reref(EEG,    []);
    EEG_TF = pop_reref(EEG_TF, []);

    % Determine the rank of data (how many linearly independent vectors there are in the data matrix)
    dataRank = sum(eig(cov(double(EEG_TF.data'))) > 1e-6);

    % Next, we epoch the EEG data. This means thatwe cut the data into small pieces around our time-locking events (the onset of the memory arrays / of the 5 circles)
    % For ERP analysis, we cut the data into epochs ranging from -300 ms to 2800 ms around memory array onset.
    [EEG, idx_to_keep] = pop_epoch(EEG, {'X'}, [-0.3, 2.8], 'newname', ['vp_', num2str(subject_id), '_epoched'], 'epochinfo', 'yes');
    EEG.trialinfo =  EEG.trialinfo(idx_to_keep, :);
    EEG = pop_rmbase(EEG, [-200, 0]);

    % For time-frequency analysis, we cut the data into epochs ranging from -800 ms to 3300 ms around memory array onset.
    [EEG_TF, idx_to_keep] = pop_epoch(EEG_TF, {'X'}, [-0.8, 3.3], 'newname', ['vp_', num2str(subject_id), '_epoched'],  'epochinfo', 'yes');
    EEG_TF.trialinfo =  EEG_TF.trialinfo(idx_to_keep, :);
    EEG_TF = pop_rmbase(EEG_TF, [-200, 0]);

    % This function detects and rejects epochs of bad data quality based on statistical properties.
    % First, we perform this trial rejection on the time-frequency dataset, as we are going to use this for independent component analysis (ICA)
    [EEG_TF, EEG_TF.rejected_epochs] = pop_autorej(EEG_TF, 'nogui', 'on');

    % Remove the rejected trials from the trialinfo array as well (so that eeg data and trial meta-data are matching)
    EEG_TF.trialinfo(EEG_TF.rejected_epochs, :) = [];

    % Run ICA on time-frequency dataset
    EEG_TF = pop_runica(EEG_TF, 'extended', 1, 'interrupt', 'on', 'PCA', dataRank);

    % Run IC-Label
    EEG_TF = iclabel(EEG_TF);

    % Find ICs that reflect (according to IC-Label) brain activity with a probability of less than 30% or that reflect eye activity with
    % a probability of more than 30%.
    EEG_TF.nobrainer = find(EEG_TF.etc.ic_classification.ICLabel.classifications(:, 1) < 0.3 | EEG_TF.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3);

    % Copy ICs to the dataset for ERP analysis. Also copy the IC-Label results and the list of ICs identified as not reflecting brain activity.
    EEG = pop_editset(EEG, 'icachansind', 'EEG_TF.icachansind', 'icaweights', 'EEG_TF.icaweights', 'icasphere', 'EEG_TF.icasphere');
    EEG.etc = EEG_TF.etc;
    EEG.nobrainer = EEG_TF.nobrainer;

    % Save the datasets with all ICs 
    pop_saveset(EEG, 'filename', ['vp_', num2str(subject_id), '_icset_erp.set'], 'filepath', PATH_ICSET, 'check', 'on');
    pop_saveset(EEG_TF, 'filename', ['vp_', num2str(subject_id), '_icset_tf.set'], 'filepath', PATH_ICSET, 'check', 'on');

    % Remove non-brain ICs from both datasets
    EEG    = pop_subcomp(EEG, EEG.nobrainer, 0);
    EEG_TF = pop_subcomp(EEG_TF, EEG_TF.nobrainer, 0);

    % Run the epoch rejection function for the ERP set as well. Remove bad epochs from the eeg data and from the trialinfo array.
    [EEG, EEG.rejected_epochs] = pop_autorej(EEG, 'nogui', 'on');
    EEG.trialinfo(EEG.rejected_epochs, :) = [];

    % Save cleaned data
    pop_saveset(EEG, 'filename', ['vp_', num2str(subject_id), '_cleaned_erp.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');
    pop_saveset(EEG_TF, 'filename', ['vp_', num2str(subject_id), '_cleaned_tf.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');

end % End subject loop


