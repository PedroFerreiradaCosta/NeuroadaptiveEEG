function [ERP_Features, ERPts] = SD_ERPFeatures(ERP_curr, TimeInfo)

% The function [ERP_Features, ERPts] = SD_ERPFeatures(ERP_curr, TimeInfo) calculates the average ERP across trials and channels. 
% It then identifies features for the Nc (latencies, amplitudes, mean amplitude). 

% Input:
% - ERP_curr: clean trials
% - TimeInfo: info with time window of interest etc

% Output: 
% - ERP_Features: structure with features for ERP components:
% Nc amplitude, Nc latency, Nc mean amplitude, confidence level,
% AUC: Nc window, confidence interval
% - ERPts: time series averaged across trials and channel of interest for later plotting

% Edits: RH; 12-09-19, adapted by ET 08-11-19; edited AG 31-03-20; last edited AG 4-01-21


%% Create ERP time series
    % select data, average across all trials and channels in 1 go
        % average and create time vector
%         ERP.ERPch = squeeze(nanmean(ERP_curr,3)); %avg across trials
    
    if strcmp(TimeInfo.ChoI, 'Fz') %select data for Fz channel
        RawChoI = squeeze(ERP_curr(4,:,:)); %select only channel 4 (Fz)
        ERPts.ERPall = squeeze(nanmean(RawChoI,2))'; 
        
    elseif strcmp(TimeInfo.ChoI, 'all') % average across all frontal channels (FC1, C1, Fz, Cz, FC2, C2)
        RawChoI = ERP_curr([2:5,7,8],:,:);
        for rr = 1:size(RawChoI,3)
            if rr == 1
                RawTrls = RawChoI(:,:,1);
            else
                RawTrls = cat(1,RawTrls,squeeze(RawChoI(:,:,rr)));
            end
            
        end  
        
        ERPts.ERPall = nanmean(RawTrls,1);

   
    else
        error('Channel of interest for Nc unrecognised in TimeInfo')
    end
    
    
    
    %% Create REF ERP time series to substract from the channel of interest 
    
% Check whether the reference channel is valid
if strcmp(TimeInfo.REFmethod,'P7')
    if ismember('P7',TimeInfo.Channels)
        REF_ind = ismember('P7',TimeInfo.Channels)==1;
        RawChoI = squeeze(ERP_curr(REF_ind,:,:));
        ERPts.ERP_REF = squeeze(nanmean(RawChoI,2))'; 
    else
        error('P7 channel as REF not present in layout')
    end
    
elseif strcmp(TimeInfo.REFmethod,'Oz')
    if ismember('Oz',TimeInfo.Channels)==1
        REF_ind = [];
            for cc = 1:8
                if strcmp(TimeInfo.Channels{cc,1}, 'Oz')
                   REF_ind = cc;
                end 
            end
        RawChoI = squeeze(ERP_curr(REF_ind,:,:));
        ERPts.ERP_REF = squeeze(nanmean(RawChoI,2))'; 

    else
        error('Oz channel not valid for REF')
    end
    
    

elseif strcmp(TimeInfo.REFmethod,'Oz&Cz')
    if ismember('Oz',TimeInfo.Channels)==1 && ismember('Cz',TimeInfo.Channels)==1
        REF_ind1 = [];
        REF_ind2 = [];
            for cc = 1:8
                if strcmp(TimeInfo.Channels{cc,1}, 'Oz')
                   REF_ind1 = cc;
                elseif strcmp(TimeInfo.Channels{cc,1}, 'Cz')
                   REF_ind2 = cc;
                end 
            end
        RawChoI1 = squeeze(ERP_curr(REF_ind1,:,:));
        RawChoI2 = squeeze(ERP_curr(REF_ind2,:,:));
        RawChoI = cat(2,RawChoI1,RawChoI2);
        ERPts.ERP_REF = squeeze(nanmean(RawChoI,2))'; 

    else
        error('Oz & Cz channels not valid for REF')
    end

elseif strcmp(TimeInfo.REFmethod,'Oz&P7')
    if ismember('Oz',TimeInfo.Channels)==1 && ismember('P7',TimeInfo.Channels)==1
        REF_ind1 = [];
        REF_ind2 = [];
            for cc = 1:8
                if strcmp(TimeInfo.Channels{cc,1}, 'Oz')
                   REF_ind1 = cc;
                elseif strcmp(TimeInfo.Channels{cc,1}, 'P7')
                   REF_ind2 = cc;
                end 
            end
        RawChoI1 = squeeze(ERP_curr(REF_ind1,:,:));
        RawChoI2 = squeeze(ERP_curr(REF_ind2,:,:));
        RawChoI = cat(2,RawChoI1,RawChoI2);
        ERPts.ERP_REF = squeeze(nanmean(RawChoI,2))'; 

    else
        error('Oz & P7 channels not valid for REF')
    end


elseif strcmp(TimeInfo.REFmethod,'all')
    RawChoI = ERP_curr;
        for rr = 1:size(RawChoI,3)
            if rr == 1
                RawTrls = RawChoI(:,:,1);
            else
                RawTrls = cat(1,RawTrls,squeeze(RawChoI(:,:,rr)));
            end
        end       
        ERPts.ERP_REF = nanmean(RawTrls,1);
   
elseif strcmp(TimeInfo.REFmethod,'none')
        ERPts.ERP_REF = zeros(1,size(ERPts.ERPall,2));
else
    error('Unrecognised REF channel(s)')
end

    % Substract the REF from the raw ERPs
    ERPts.ERPraw = ERPts.ERPall;
    ERPts.ERPall = ERPts.ERPall - ERPts.ERP_REF;
    
    % Add time and choi to ERPts structure    
    ERPts.Time_ms = TimeInfo.Time*1000;
    ERPts.ChoI = TimeInfo.ChoI;

   
%% Nc amplitude, latency and mean amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find peak Nc
    PeakNc = find(ERPts.ERPall(1,TimeInfo.ToI_Nc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)) == min(ERPts.ERPall(1,TimeInfo.ToI_Nc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)))); 
    Nc_Ampl = ERPts.ERPall(TimeInfo.ToI_Nc_samps(1,1)+PeakNc-1);
    Nc_MAmpl = mean((ERPts.ERPall(1,TimeInfo.ToI_Nc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2))),2); % 2 for 2nd dimension ( i.e. columns/time points), 1 for rows program: sum up all these amplitudes embraced in these brackets     
   

    ERP_Features.Nc_Lat = ERPts.Time_ms(TimeInfo.ToI_Nc_samps(1,1) + PeakNc-1);
    ERP_Features.Nc_Ampl = Nc_Ampl;
    ERP_Features.Nc_MAmpl = Nc_MAmpl;
    
    
%% Calculate confidence levels 
    % Select baseline interval for ERPall
        BaselineERP = ERPts.ERPall(1,TimeInfo.ToCi_samps(1,1):TimeInfo.ToCi_samps(1,2));
    % Calculate confidence level (standard deviation during baseline)
        ERP_Features.Conf_stdev = std(BaselineERP); %not relevant?   
        
        
        
        
%% AUC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create abs values
    AbsERP = abs(ERPts.ERPall);
  
    % Nc timewindow
    ERP_Features.AUC_Nc_window = trapz(TimeInfo.Time(1,TimeInfo.ToI_Nc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)),AbsERP(1,TimeInfo.ToI_Nc_samps(1,1):TimeInfo.ToI_Nc_samps(1,2)));
     
    % confidence interval
    ERP_Features.AUC_ConInt = trapz(TimeInfo.Time(1,TimeInfo.ToCi_samps(1,1):TimeInfo.ToCi_samps(1,2)),AbsERP(1,TimeInfo.ToCi_samps(1,1):TimeInfo.ToCi_samps(1,2))); %not relevant?   
        
        
        
        
        
        
        end        

 
        

