function [FigERPs] = SD_ERPts_Figure(FigERPs, ERPs, ERP_Features, Ntrials, Blck, TimeInfo, face_N)

% The function [FigERPs] = SD_ERPts_Figure(FigERPs, ERPs, ERP_Features, Blck, TimeInfo)
% plots the ERP time series and markers at the ERP components of interest.

% Input:
% - FigERPs: figure to replot new ERP in
% - ERPs; time series averaged across trials and channels
% - ERP_Features: structure with features ERP components
% - Ntrials: number of trials per block
% - Blck: number of the current block
% - TimeInfo: info with time window of interest etc
% - face_N: marker (face 1 or face 2)

% Output: 
% - FigERPs: figure with new ERP 

% RH; 12-09-19, updated by AG 02/12/2019
   
%% Create strings for the title    
   
% create string for Nc measures
    formatSpec = 'Nc: %.2f uV, %.2f ms, %.2f uV, Quality (%.2f)';
    Rep_Nc = sprintf(formatSpec, ERP_Features.Nc_Ampl, ERP_Features.Nc_Lat, ERP_Features.Nc_MAmpl, ERP_Features.Conf_stdev);
       
%% plot ERP with markers at identified peaks

% get the figure and prep
    figure(FigERPs)
    clf %clear current figure window
    set(FigERPs, 'visible','off') %make invisible and plot 
    
% plot the data and peak    
    plot(ERPs.Time_ms,ERPs.ERPall,'b')
    hold on
    
    % plot Nc peak
    plot(ERP_Features.Nc_Lat,ERP_Features.Nc_Ampl,'Marker','.','MarkerSize',20) 
    ylim([-25 25])
    
    
    
% plot the labels and axes through the origin  
    xlabel('Time (ms)'); ylabel('Amplitude (uV)')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    
    
    
% plot the title with info    
    Rep_title = ['ERPs Face ' face_N ' (' num2str(Ntrials) ' trials in block ' num2str(Blck) ') for ' TimeInfo.ChoI];
    title({Rep_title; ' '; Rep_Nc}, 'Fontsize',12);

    % set figure to visible again    
    set(FigERPs, 'visible','on') %make plot visible
    
    

 
end
