function [FigChQual, Incl] = SD_Channel_Feedback(Ch_curr, TimeInfo, FigChQual, Ntrials, Blck, N_goodQualERP_aim, perc_good)

% The function [FigChQual] = SD_Channel_Feedback(Ch_curr, TimeInfo, FigChQual)
% plots the feedback for each channel, number of trials surviving artefact 
% rejection AR threshold, range, and flat signals.

% Input:
% - Ch_curr; structure with channel AR info (set to NaN across time for that trial in ERP_curr) 
%   - Excl: excl = 1, incl = 0 in trials, ch x trl
%   - SumARs: info on ARmethod for each channel, ch x ARmethod (thr, ran, flat) x trl
% - TimeInfo: info with time window of interest etc
% - FigChQual: figure to plot feedback in

% Output: 
% - FigChQual: figure with feedback for current block

% Adapted from RH GxE_Channel_Feedback 13-2-20; last update AG 11-2-21

%% Calculate values for bar plot
% for each channel; N trials included - N out by thr - range - flat 
% ch x Ns


    Channels = TimeInfo.Channels;
    N_trls_curr = size(Ch_curr.Excl,2); %remember to change all Ch_curr.Excl into Ch_curr{Blck}.Excl if outside realtime loop
    Incl = zeros(size(Channels,1),4);
    
    if strcmp(TimeInfo.Protocol, 'OnlineEEG_wholescalp_wLS')
        % include N trials incl for LS
        Incl(1,:) = repmat(N_trls_curr, [1,4]);
        % include N trials incl for other channels
        % all
        Incl(2:8,1) = N_trls_curr - sum(Ch_curr.Excl,2);
        % threshold
        Incl(2:8,2) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,1,:)),2);
        % range
        Incl(2:8,3) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,2,:)),2);
        % flat
        Incl(2:8,4) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,3,:)),2);
    else
        % all
        Incl(1:8,1) = N_trls_curr - sum(Ch_curr.Excl,2);
        % threshold
        Incl(1:8,2) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,1,:)),2);
        % range
        Incl(1:8,3) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,2,:)),2);
        % flat
        Incl(1:8,4) = N_trls_curr - sum(squeeze(Ch_curr.SumARs(:,3,:)),2);
    end

    disp('Wait for channel feedback')
        figure(FigChQual)
        clf 
        set(FigChQual, 'visible','off')  
        % create good/bad areas
        rectangle('Position',[0 0 9 N_goodQualERP_aim],'FaceColor',[0.6350 0.0780 0.1840 .2], 'EdgeColor','none');
        hold on
        rectangle('Position',[0 N_goodQualERP_aim 9 (Ntrials-N_goodQualERP_aim)],'FaceColor',[0.4660 0.6740 0.1880 .2], 'EdgeColor','none');
        % plot bars
        b = bar(Incl);
        % change the colors of the bars
        b(1).FaceColor = 'k'; % total
        b(2).FaceColor = [0 0.4470 0.7410]; % thr
        b(3).FaceColor = [0.4940 0.1840 0.5560]; % range
        b(4).FaceColor = [0.3010 0.7450 0.9330]; % flat
        % add labels
        xticks([1 2 3 4 5 6 7 8])
        xticklabels(Channels')
        xlabel('Channels','FontSize',14)
        ylabel('Number of trials included','FontSize',14)
        legend({'Total','Within thresholds','Within range','No flat signal'},'Orientation','horizontal','Location','southoutside')
        perc_good_str = string(perc_good)
        title({'Trials included per channel'; strcat('Block ',num2str(Blck), ' | % of total time series (6 channels X 12 trials) included in ERP: ', perc_good_str)}, 'Fontsize',14); 
        set(FigChQual, 'visible','on') %make plot visible  
    disp('Plotting done')

end