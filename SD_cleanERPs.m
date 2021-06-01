function [ERP_curr, Ch_curr] = SD_cleanERPs(eeg, eeg_t, t_mark, TimeInfo)

% The function [ERP_curr, Ch_curr] = SD_cleanERPs(eeg, eeg_t, t_mark, TimeInfo, ARInfo) 
% preprocesses the raw EEG LSL data for later ERPs averaging.

% - filter: band pass .1 - 20 Hz, 499th order, Hanning window
% - segmenting: segmenting trials based on TimeInfo
% - baseline correction: based on TimeINfo
% - artefact detection: threshold, range, flat for time of interest for
% component
% Across all channels

% Input:
% - eeg; raw eeg data from online enobio
% - eeg_t; timesamples for each eeg sample
% - t_mark; timesamples for each eeg marker
% - TimeInfo: structure with info on samples, output from SD_Timings script

% Output: 
% - ERP_curr; preprocessed clean trials, ch x time x trls (excl data set to
% NaN)
% - Ch_curr; structure with channel AR info (set to NaN across time for that trial in ERP_curr) 
% - Excl: excl = 1, incl = 0 in trials, ch x trl
% - SumARs: info on ARmethod for each channel, ch x ARmethod (thr, ran, flat) x trl


%% Filtering on all trials

% prep for data filtering

    EEGraw2 = eeg; % in uV data
    % demean
    EEGmn = mean(EEGraw2,2);
    EEGrawdm = EEGraw2 - repmat(EEGmn,[1,size(EEGraw2,2)]);
    % transpose for filtering (along first >1 dimension)
    EEGrawdmprep = EEGrawdm'; 
    

% Finite response filter
% prep filter
    BPfilter = fir1(499,[1/250 20/250],'bandpass',hann(500)); %bp from 1 to 20 Hz, Fs = 500, Fn = 250 Hz
    
% apply filter to the prep data
    FiltEEG = filtfilt(BPfilter,1,EEGrawdmprep);
% transpose data chxtime
    FiltEEG = FiltEEG';
       


%% Segmenting and baseline correction on all channels and trials


NTrls = length(t_mark); %Number of trials/markers in the block
TrialData = zeros(size(FiltEEG,1),(sum(TimeInfo.TSeg_samps,2)+1),NTrls);


Begs_samples = zeros(NTrls,1);
Ends_samples = zeros(NTrls,1);

for bb = 1:NTrls
    A = find((eeg_t - t_mark(1,bb))>= 0, 1);
    BEG = A - TimeInfo.TSeg_samps(1,1);
    END = A + TimeInfo.TSeg_samps(1,2);
    Begs_samples(bb,1) = BEG;
    Ends_samples(bb,1) = END;
    clear BEG END
end
clear bb

% check end samples
    Valid_beg_ind = (Begs_samples >= 1); %check begin earlier than 1
    Valid_end_ind = (Ends_samples <= size(FiltEEG,2)); % check end later than recorded eeg samples
    
    Valid_both = Valid_beg_ind+Valid_end_ind; % sum beg and end valid, should be 2 if both valid
    Ind_validboth_ind = (Valid_both == 2); % find out  which trials to select if valid
    
    Beg_validall = Begs_samples(Ind_validboth_ind,1); % select valid begin samples
    End_validall = Ends_samples(Ind_validboth_ind,1); % select valid end samples
    
% segment data
for cc = 1:size(Beg_validall,1) % segment trials based on valid begin and end samples
    BEG = Beg_validall(cc,1); % find current begin sample
    END = End_validall(cc,1); % find current end sample
    TrialData(:,:,cc) = FiltEEG(:,BEG:END); % select data and save into trial data
end
clear cc

    


% for bb = 1:NTrls
%     
%     A = find((eeg_t - t_mark(1,bb))>= 0, 1);
%     if isempty(A)
%         % handle no data returned => first t_mark is earlier than t_eeg! or
%         % last t_mark later than t_eeg (last trial)
%     end
%     
%     BEG = A - TimeInfo.TSeg_samps(1,1);
%     END = A + TimeInfo.TSeg_samps(1,2);
%     
%     if END > size(FiltEEG,2)
%         % handle no data returned => last t_mark is later than t_eeg! only
%         % error for last trial in block
%         %if >, look at t for previous one, look if that one is smaller, if
%         %not look at t for next previous one
%     end
%     TrialData(:,:,bb) = FiltEEG(:,BEG:END); 
% end

% clear bb

%% Baseline correction 

% calculate mean across baseline time 
    BlAvg = mean(TrialData(:,TimeInfo.TBc_samps(1,1):TimeInfo.TBc_samps(1,2),:),2); % ch x 1 x NTrls
% create matrix with baseline values to subtract
    Bl = repmat(BlAvg,[1,(sum(TimeInfo.TSeg_samps,2)+1),1]);
% apply baseline correction
    BlcorrSeg = TrialData - Bl;

%% Artefact detection and average across channels for each trial separately
% only take out channel when artefacts occur during the timewindow of interest ToI

% pre-allocate data
[Nch, Ntime, ~] = size(BlcorrSeg);
ERP_curr = zeros(Nch, Ntime, NTrls); %ch x time x trl
Ch_curr.Excl = zeros(Nch, NTrls); %ch x trl
Ch_curr.SumARs = zeros(Nch, 3, NTrls); %ch x ARmethod (thr, ran, flat) x trl

% identify and remove artefacts for each trial
for rr = 1:NTrls

    TrlCur = BlcorrSeg(:,:,rr);

    % channels exceeding thresholds
    ARthr = any(TrlCur(:,TimeInfo.TBc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)) < TimeInfo.ARInfo.thrs(1,1),2)| any(TrlCur(:,TimeInfo.TBc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)) > TimeInfo.ARInfo.thrs(1,2),2); 
    % channels exceeding range
    ARran = max(TrlCur(:,TimeInfo.TBc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)),[],2) - min(TrlCur(:,TimeInfo.TBc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)),[],2) >= TimeInfo.ARInfo.range;  
    % channels with flat signals
    ARflat = all(abs(TrlCur(:,TimeInfo.TBc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2))) < .0001,2);

    % summarise ARs 
    SumARs = [ARthr, ARran, ARflat];

    % create NaNs for the channel with artefacts
    TbExcl = sum(SumARs,2);
    TrlCur(TbExcl>=1,:) = NaN;

    % save everything in the variables
    ERP_curr(:,:,rr) = TrlCur; %put trialdata in variable
    Ch_curr.Excl(:,rr) = TbExcl; %put channel incl 0/ excl 1 in variable
    Ch_curr.SumARs(:,:,rr) = SumARs; %summary of ARs: ARthr ARran ARflat

end

end