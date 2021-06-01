function [TimeInfo] = SD_Timings(TimeInfo)

% This function derives the indices and timing for later segmeting of the EEG data for the ERPs during Faces Space Definition with Nc (SD_Nc_faces).

% Input:
% - TimeInfo: structure with timings in seconds for:
% segmenting, baseline correction, confidence interval,
% and time of interest (in seconds) for Nc
% and EEG protocol

% Output: 
% - TimeInfo: structure with timing info and samples for segmenting,
% baseline, window of interest, window for confidence levels, and list of EEG channels for the
% current protocol. 

% last edited AG 02/12/2019 based on RH's scripts (last edit 18-11-19):

%% Check inputs

    if ~isequal(size(TimeInfo.TSeg,2),2) 
        error('Info missing; check times for segmenting')
    end

    
    if ~isequal(size(TimeInfo.TBc,2),2) 
        error('Info missing; check times for baseline correction')
    end

    
    if ~isequal(size(TimeInfo.ToCi,2),2) 
        error('Info missing; check times for confidence interval')
    end

    
    if ~isequal(size(TimeInfo.ToI_Nc,2),2)
        error('Info missing; check times for window of interest for Nc')
    end

    
    if isempty(TimeInfo.Fs)
        error('Info missing; check sampling rate')
    end

    
    if isempty(TimeInfo.Protocol)
       error('Info missing; specify protocol')
    end

    
    if isempty(TimeInfo.ARInfo.thrs)
       error('Info missing; specify threshold for AR')
    end

    
    if isempty(TimeInfo.ARInfo.range)
       error('Info missing; specify range for AR')
    end

    
    if isempty(TimeInfo.REFmethod)
       error('Info missing; reference method is missing')
    end

%% Derive indices

TimeInfo.Time = TimeInfo.TSeg(1,1):1/TimeInfo.Fs:TimeInfo.TSeg(1,2);
TimeInfo.Marker_t0 = find(TimeInfo.Time == 0);

% Segment samples around marker t=0
    TimeInfo.TSeg_samps(1,1) = find(TimeInfo.Time == 0)-1;
    TimeInfo.TSeg_samps(1,2) = size(TimeInfo.Time,2) - find(TimeInfo.Time == 0);

% Baseline correction samples within segmented time
    if  TimeInfo.Time(1,1) > TimeInfo.TBc(1,1) || TimeInfo.TBc(1,2) > TimeInfo.Time(1,end) % check if specified baseline is valid
        error('Invalid time for baseline correction (earlier or later sample than segmented data).')
    else
        TimeInfo.TBc_samps(1,1) = find(abs(TimeInfo.Time - TimeInfo.TBc(1,1)) < .001);
        TimeInfo.TBc_samps(1,2) = find(abs(TimeInfo.Time - TimeInfo.TBc(1,2)) < .001);
    end

% Time window of interest in samples within segmented time
    if  TimeInfo.Time(1,1) > TimeInfo.ToCi(1,1) || TimeInfo.ToCi(1,2) > TimeInfo.Time(1,end) % check if specified window is valid
        error('Invalid time for window of confidence (earlier or later sample than segmented data).')
    else
        TimeInfo.ToCi_samps(1,1) = find(abs(TimeInfo.Time - TimeInfo.ToCi(1,1)) < .001);
        TimeInfo.ToCi_samps(1,2) = find(abs(TimeInfo.Time - TimeInfo.ToCi(1,2)) < .001);
    end

%% ERP components; windows of interest


% Time window of interest in samples within segmented time: Nc
    if  TimeInfo.Time(1,1) > TimeInfo.ToI_Nc(1,1) || TimeInfo.ToI_Nc(1,2) > TimeInfo.Time(1,end) % check if specified window is valid
        error('Invalid time for Nc window of interest (earlier or later sample than segmented data).')
    else
        TimeInfo.ToI_Nc_samps(1,1) = find(abs(TimeInfo.Time - TimeInfo.ToI_Nc(1,1)) < .001);
        TimeInfo.ToI_Nc_samps(1,2) = find(abs(TimeInfo.Time - TimeInfo.ToI_Nc(1,2)) < .001);
    end


%% Channels and layout checks

% Channel order
    if strcmp(TimeInfo.Protocol,'OnlineEEG_wholescalp')
        TimeInfo.Channels = {'P7';'FC1';'C1';'Fz'; 'Cz'; 'Oz';'FC2';'C2'};
    elseif strcmp(TimeInfo.Protocol,'OnlineEEG_wholescalp_wLS')
        TimeInfo.Channels = {'LS';'FC1';'C1';'Fz'; 'Cz'; 'Oz';'FC2';'C2'};
    else
        error('Protocol not recognized')
    end



% Check whether channel/s of interest is valid
    if strcmp(TimeInfo.ChoI,'Fz')
        if ~ismember(TimeInfo.Channels, TimeInfo.ChoI)
            error('Fz channel of interest not present in layout')
        end
    elseif strcmp(TimeInfo.ChoI,'all')
        if ~isequal(sum(ismember(TimeInfo.Channels, {'FC1';'C1';'Fz';'Cz';'FC2';'C2'})),6)
            error('Frontal channels of interest not present in layout')
        end    
    end



end
