%% Nc_gazeXemotion: Nc for faces with direct vs averted gaze and neutral to smiling face

% Adapted from RH's script for GANxEEG experiment for the presentation of faces and
% analyses of N170 peaks. AG last updated: 11 apr 2021



% Set up for the task  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 1) open NIC 
% 2) connect to NIC box via wifi outside of Matlab
% 3) enable synchronizer
% 4) ensure only 1 stimulus is in the stimulus folder

% then run script:
clear all
close all
addpath(genpath('/Users/braintools/Documents/MATLAB/TaskEngine2'));     
addpath('/Users/braintools/Documents/MATLAB/lm_tools');
addpath('/Users/braintools/Desktop/BONDS_SD/');
addpath('/Users/braintools/Desktop/BONDS_SD/stimulus');
addpath('/Users/braintools/Desktop/BONDS_SD/data');


% Define parameters for task

NumScreen = 1; %number of monitors; change to 2 in 202, change to 1 if only 1 screen attached

% Presentation timings: 
% 300-500 fixation cross, 500 duration stim, get rid of blank in between
    dur_fix_min = 0.300;
    dur_fix_max = 0.500;
    dur_stim = 0.500;
    
% Numbers of trials and blocks
    Ntrials = 12; %number of trials within block
    NBlcks = 15; %number of blocks

% Numbers of good ERPs
    N_goodQualERP_aim = 5;    
    
% EEG processing parameters    
    % EEG processing parameters    
    TimeInfo.TSeg = [-.1 .8]; %segmenting in seconds
    TimeInfo.TBc = [-.1 0]; %baseline correction in seconds
    TimeInfo.ToCi = [-.1 0]; %timewindow confidence interval
   
    TimeInfo.ToI_Nc = [.3 .8]; %timewindow of interest for Nc; apply SD_findTW function below to adapt timewindow individually per block
    
    TimeInfo.Fs = 500; % sampling rate in Hz
    TimeInfo.Protocol = 'OnlineEEG_wholescalp'; % 'OnlineEEG_wholescalp' 'OnlineEEG_wholescalp_wLS'
    TimeInfo.ChoI = 'all'; %channel of interest; 'Fz' or 'all'
    TimeInfo.REFmethod = 'Oz&P7'; %'all','Oz', 'Cz','Oz&Cz', 'none', 'P7', 'Oz&P7'
   
 
% EEG preprocessing: AR parameters
    TimeInfo.ARInfo.thrs = [-200, 200];
    TimeInfo.ARInfo.range = 400;

% EEG buffer time (in sec, for a) EEG marker buffer, and b) padding for
% filtering)
    BufferTime = 3;

% Name ppt and create new folder
    prompt = 'Participant ID: ';
    Subj_ID = input(prompt,'s');
    if isempty(Subj_ID)
        Subj_ID = 'XXX';
    end
    mess = ['Starting session for: ', Subj_ID];
    disp(mess)
    
    Name_pptfolder = strcat('/Users/braintools/Desktop/BONDS_SD/data/',Subj_ID);
    if isfolder(Name_pptfolder)==0
        mkdir(Name_pptfolder)
    end
    
 
% output info  
    Outputmeasure = 'Nc_MAmpl'; % options: 'Ampl'; 'Lat'; Nc_MAmpl'
    
% create empty output files    
    dlmwrite('/Users/braintools/Desktop/BONDS_SD/Output.txt',[])   
     
    disp('Ready for next section')
    
% %%%%%%%% check monitor number %%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    prompt = 'Start Python, then press c to continue: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'c';
    end

% Set-up parameters for EEG

% set counters to 0 and create empty variables
    Blck = 1;
    ClnTrials = [];
    
% Get the timing info
    [TimeInfo] = SD_Timings(TimeInfo);
    
% open figures for later plotting
    Fig_Nc = figure;
    set(gcf,'position',[896, 572, 785, 383],'MenuBar','none')
    plot(0,0,'k','Marker','.','MarkerSize',20)
    
    FigChQual = figure; % Channel feedback
    set(gcf,'position',[1112 385 781 330],'MenuBar','none')
    plot(0,0,'k','Marker','.','MarkerSize',20)
    
% Name for data file
     NameMat = strcat(Name_pptfolder, '/Faces_eeg_output.mat');    

% set up task engine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pres = tePresenter;
    pres.SkipSyncTests = true;
    pres.MonitorNumber = NumScreen;
    pres.SetMonitorDiagonal(24, 16, 9, 'inches') % diagonal in '', ratio width: height
    pres.PreviewScale = .15; %size of preview screen
    pres.PreviewPositionPreset = 'bottomright';  %location of preview screen
    pres.OpenWindow
    pres.BackColour = [120, 120, 120]; %background color: 0 0 0 = white, 255 255 255 = black
%     pres.LightPatchEnabled = true;  % set up LS
%     pres.LightPatchSize = [50,50]; %size in pixels
    
    
    % set PTB font size for fixation cross
    Screen('TextSize', pres.WindowPtr, 70); 
    
    % set up event type for EEG
    pres.EventRelays('eeg') = teEventRelay_Enobio;
    
    % load attention getter images
    pres.LoadStim('/Users/braintools/Desktop/BONDS_SD/SD_attention_getters')
    

% Load stimuli into presenter from the stimulus folder
    gs = ganStimGxE('/Users/braintools/Desktop/BONDS_SD/stimulus', 1, pres);
    numLoaded = gs.LoadNewFilesIntoPresenter;
    if numLoaded > 1
        error('Loaded more than one image from the folder - was only expecting 1.')
    end
%  %%
% % connect to LSL (= "lab screen layer", this is a default code provided by the company) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load lsl library
    lib = lsl_loadlib;
    
    outletName = 'LSLOutletStreamName-EEG';
    result = {};
    while isempty(result)
        result = lsl_resolve_byprop(lib, 'type', 'EEG'); 
    end

    index=-1;
    for r=1:length(result)
        if (strcmp(name(result{r}),outletName)==1)
            index=r;
            disp('NIC stream available')
        end
    end

    if (index == -1)
        disp('Error: NIC stream not available \n');
        return;
    end

    disp('Connecting to NIC stream...');
    inlet_eeg = lsl_inlet(result{index});
    
    
% NIC markers
    outletName = 'LSLOutletStreamName-Markers';
    result = {};
    while isempty(result)
          result = lsl_resolve_byprop(lib,'type','Markers'); 
    end

    index=-1;
    for r=1:length(result)
        if (strcmp(name(result{r}),outletName)==1)
            index=r;
            disp('NIC stream available')
        end
    end

    if (index == -1)
        error('Error: NIC stream not available \n');
    end

    disp('Connecting to NIC stream...');
    inlet_marker = lsl_inlet(result{index});   

   disp('Press SPACE while child looks at screen');


% Start the experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disable keyboard
ListenChar(2)

while Blck <= NBlcks 

    % disable keyboard
 %       ListenChar(2)

%%%% this is just to replace python script

% % delete cached (.baked) files from last block
         delete('/Users/braintools/Desktop/BONDS_SD/stimulus/*.baked.mat')
 
        
%%%% if python ready, then run from here
    
% Set up EEG for this block
    % set all parameters to 0
        countPres = 0;
        inlet_marker.pull_sample(0); %empty marker buffer
        inlet_eeg.pull_chunk(); %empty eeg buffer
        
    % wait BufferTime
    pause(BufferTime) 

    
% Presentation of the stimuli for Ntrials


    for ii = 1:Ntrials
        
        if isequal(ii,1)
            dur_fix = 1;
        else
        % calculate jitter fixation period for this trial
            dur_fix_range = dur_fix_max - dur_fix_min;
            dur_fix = dur_fix_min + (rand * dur_fix_range);
        end
        
     

        % fixation 500ms
        DrawFormattedText(pres.WindowPtr, '+', 'center', 'center', [255, 255, 255]);
        pres.RefreshDisplay;
        WaitSecs(dur_fix);
        
        % present stimulus
        checkKeyboard(pres);
            
        % darw stim and wait 500ms
        pres.DrawStim('ganstim01');
%         pres.LightPatchOn;
        onset = pres.RefreshDisplay;
        pres.SendEvent(1, onset);
        WaitSecs(dur_stim); %presentation time of stimulus in sec
       
        % add 1 to counter for trials
        countPres = countPres+1;
    
    end
    
        % re-enabled keyboard
        ListenChar(0)%if keyboard gets stuck, press Control c

        % present attention getters between blocks 
        breakimg = pres.Stim.LookupRandom('Key', 'fix_img_baby*'); %  to do: load attention getters into folder and name them like this
        pres.DrawStim(breakimg);
        pres.RefreshDisplay;

        % wait BufferTime 
        pause(BufferTime)     
        
        
        
% EEG processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get markers and time samples
        [mark, t_mark] = inlet_marker.pull_chunk(); 
    % get EEG data and time samples
        [eeg, eeg_t] = inlet_eeg.pull_chunk();
            
    % preprocess trials
        [Pre_trials, ChExcl] = SD_cleanERPs(eeg, eeg_t, t_mark, TimeInfo);

    % calculate ERP across channels of interest
        [ERP_Features, ERPs] = SD_ERPFeatures(Pre_trials, TimeInfo);
        
       
%             %if at least half of the channels give enough trials for ERPs, get and save output for python...    
%     
%         if strcmp(TimeInfo.ChoI,'Fz')
%            ChIndex=find(strcmp(TimeInfo.Channels,'Fz'));
%         elseif strcmp(TimeInfo.ChoI,'all')
%            ChoIonly = {'FC1' 'C1' 'Fz' 'Cz' 'FC2' 'C2'};
%            for a = 1:length(ChoIonly)
%                 Ch_match = strcmp(TimeInfo.Channels,ChoIonly(a));
%                 ChIndex(a) = find(Ch_match);
%            end 
%         end
%            
%           
%         GoodCh = sum(Incl(ChIndex,1) > N_goodQualERP_aim);
       
 % check timeseries for 12 trials X 6 channels (WITHOUT references)

            if strcmp(TimeInfo.ChoI,'Fz')
               ChIndex=find(strcmp(TimeInfo.Channels,'Fz'));
            elseif strcmp(TimeInfo.ChoI,'all')
               ChoIonly = {'FC1' 'C1' 'Fz' 'Cz' 'FC2' 'C2'};
               for a = 1:length(ChoIonly)
                    Ch_match = strcmp(TimeInfo.Channels,ChoIonly(a));
                    ChIndex(a) = find(Ch_match);
               end 
            end 

        % drop reference rows manually (Ch_Excl_8 = incl ref; Ch_Excl_6 = without reference electrodes)
        ChExcl_8 = ChExcl       
        ChExcl.Excl(6,:) = [] 
        ChExcl.Excl(1,:) = []
        ChExcl_6 = ChExcl
        ChExcl = ChExcl_8
        count = sum(ChExcl_6.Excl(:) == 1)
        total = size(ChExcl_6.Excl(:))
        perc_bad = (count/total(1))*100
        perc_good = 100 - perc_bad


%     % check timeseries for 12 trials X 8 channels (INCL references)  
%
%      count = sum(ChExcl.Excl(:) >= 1)
%      total = size(ChExcl.Excl(:))
%      perc_bad = (count/total(1))*100
%      perc_good = 100 - perc_bad

   
[FigChQual,Incl] = SD_Channel_Feedback_perc(ChExcl_8, TimeInfo, FigChQual, Ntrials, Blck, N_goodQualERP_aim, perc_good);


%% check data qual and decide whether to repeat

% if GoodCh>length(ChIndex)/2 %at least half of the channels have enough trials
    prompt = 'Check data quality - continue (c), repeat block (r) or switch to offline (o)?';
    str1 = input(prompt,'s');

    if contains(str1,'c')
        
        
        % Save outputmeasure to txt files  
            [OutputVal] = SD_ERPoutputVal(ERP_Features, Outputmeasure);
            disp('Txt files written')    
  
        % python can start looking and generating now while current script plots output and saves the data in variables


        % Report and bookkeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot final ERP with identified features
            Ntrials = size(Pre_trials,3);
            [Fig_Nc] = SD_ERPts_Figure(Fig_Nc, ERPs, ERP_Features, Ntrials, Blck, TimeInfo, ' ');

        % EEGraw data into structure for this block
            EEGraw = struct();
            EEGraw.eeg = eeg;
            EEGraw.eeg_t = eeg_t;
            EEGraw.mark = mark;
            EEGraw.t_mark = t_mark;    

        % append clean trials and latencies
            if isequal(Blck,1) % for the first block
                % clean eeg data
                    ClnTrials = Pre_trials;
                % Nc
                    Output.Nc_Lat = ERP_Features.Nc_Lat;
                    Output.Nc_Ampl = ERP_Features.Nc_Ampl;
                    Output.Nc_MAmpl = ERP_Features.Nc_MAmpl;
                % For Bayesian optimasation online
                    Output.EEGoutputVal = OutputVal;

                % ERP time series
                    ERPtimeseries{1,1} = ERPs;
                % Channel and eeg data
                    Ch_Excl{1,1} = ChExcl;
                    EEGdata_raw{1,1} = EEGraw;

            else % concatenate with other blocks
                % clean eeg data
                    ClnTrials = cat(4,ClnTrials,Pre_trials); %Ch x Time x Trl x Blck            
                % Nc
                    Output.Nc_Lat = cat(1,Output.Nc_Lat,ERP_Features.Nc_Lat);
                    Output.Nc_Ampl = cat(1,Output.Nc_Ampl,ERP_Features.Nc_Ampl);
                    Output.Nc_MAmpl = cat(1,Output.Nc_MAmpl,ERP_Features.Nc_MAmpl);
                % For Bayesian optimasation online
                    Output.EEGoutputVal = cat(1,Output.EEGoutputVal,OutputVal);
                % ERP time series    
                    ERPtimeseries = cat(1,ERPtimeseries,ERPs); %Blck
                % Channel and eeg data
                    Ch_Excl = cat(1,Ch_Excl,ChExcl); %Blck
                    EEGdata_raw = cat(1,EEGdata_raw,EEGraw); %Blck
            end

            clear Pre_trials ChExcl 
            clear ii countPres eeg eeg_t mark t_mark
            clear ERP_Features

            NameFig = strcat(Name_pptfolder,'/ERPs_Blck_',num2str(Blck),'.png');
            saveas(Fig_Nc, NameFig)

            NameMatBlock = strcat(Name_pptfolder, '/Faces_eeg_output_blck',num2str(Blck),'.mat');    
            save(NameMatBlock)

            % Await new stimuli from GAN and load them in  
            
            if Blck < NBlcks        

                fprintf('Waiting for new stimuli...\n')
                gs.FolderHasChanged;
                while gs.NumNewFiles < 1
                    str1 = [];
                    if isempty(str1)
                        prompt = 'Check Python - Continue session? y/n';
                        str1 = input(prompt,'s');
                         if contains(str1,'n')

                            Faces = Output;
                            save(NameMat,'ClnTrials','Faces','Ch_Excl', 'EEGdata_raw','TimeInfo','Output')
                            save(strcat(Name_pptfolder,'/EEGdata_raw.mat'), 'EEGdata_raw')

                            close all
                            clear all

                            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                            disp('End of the session (py)')
                            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')     
                            return
                            
                         end 
                    end
                  gs.FolderHasChanged;
                  WaitSecs(0.01);
                end 
                
                WaitSecs(8);
                gs.RemoveOldFilesFromPresenter
                gs.LoadNewFilesIntoPresenter

                endblock = sprintf('End of block %2d \n', Blck);
                disp(endblock);

                % save variables into .mat file at end of the session
                Faces = Output;
                save(NameMat,'ClnTrials','Faces','Ch_Excl', 'EEGdata_raw','TimeInfo','Output')
                save(strcat(Name_pptfolder,'/EEGdata_raw.mat'), 'EEGdata_raw')

                
                Blck = Blck + 1;           
            end
                         
        elseif contains(str1,'r')
                    rep_text = sprintf('Repeating block %2d',Blck);
                    disp(rep_text);
                    Blck = Blck;%if not enough ts with good data      
            
    else 
            
        OutputVal = 999;  %save 999 in Output file for py
                    dlmwrite('/Users/braintools/Desktop/BONDS_SD/Output.txt',OutputVal);

                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp('Continuing without BO, showing predefined images')
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

                    dlmwrite('/Users/braintools/Desktop/BONDS_SD/data/LastBlockOnline.txt',Blck );
 %                 return
  
 addpath('/Users/braintools/Desktop/BONDS_SD/9points/images');
 pres.LoadStim('/Users/braintools/Desktop/BONDS_SD/9points/images')

                    i = 0
                    
                    while Blck <= NBlcks                                                
                        
                    % Set up EEG for this block
                        % set all parameters to 0
                        countPres = 0;
                        inlet_marker.pull_sample(0); %empty marker buffer
                        inlet_eeg.pull_chunk(); %empty eeg buffer

                    % wait BufferTime
                        pause(BufferTime) 

                    % Presentation of the stimuli for Ntrials
                        for ii = 1:Ntrials

                            if isequal(ii,1)
                                dur_fix = 1;
                            else
                            % calculate jitter fixation period for this trial
                                dur_fix_range = dur_fix_max - dur_fix_min;
                                dur_fix = dur_fix_min + (rand * dur_fix_range);
                            end

                    % fixation 500ms
                            DrawFormattedText(pres.WindowPtr, '+', 'center', 'center', [255, 255, 255]);
                            pres.RefreshDisplay;
                            WaitSecs(dur_fix);

                    % present stimulus
                            checkKeyboard(pres);

                     % draw stim and wait 500ms pres.DrawStim('ganstim01');

                            if i == 0
                                pres.DrawStim('tmp00.png');
                            elseif i == 1
                                pres.DrawStim('tmp15.png');
                            elseif i == 2
                                pres.DrawStim('tmp12.png');
                            elseif i == 3
                                pres.DrawStim('tmp03.png');
                            elseif i == 4
                                pres.DrawStim('tmp04.png');
                            elseif i == 5
                                pres.DrawStim('tmp13.png');
                            elseif i == 6
                                pres.DrawStim('tmp11.png');
                            elseif i == 7
                                pres.DrawStim('tmp02.png');
                            elseif i == 8
                                pres.DrawStim('tmp08.png');
                            elseif i == 9
                                pres.DrawStim('tmp14.png');
                            elseif i == 10
                                pres.DrawStim('tmp07.png');
                         elseif i == 11
                                pres.DrawStim('tmp01.png');
                            elseif i == 12
                                pres.DrawStim('tmp09.png');
                            elseif i == 13
                                pres.DrawStim('tmp06.png');
                            elseif i == 14
                                pres.DrawStim('tmp10.png');
                            else 
                                pres.DrawStim('tmp05.png');
                            end

                            onset = pres.RefreshDisplay
                            pres.SendEvent(1, onset);
                            WaitSecs(dur_stim); %presentation time of stimulus in sec

                    % add 1 to counter for trials
                            countPres = countPres+1;

                      end

                    % present attention getters between blocks 
                        breakimg = pres.Stim.LookupRandom('Key', 'fix_img_baby*'); %  to do: load attention getters into folder and name them like this
                        pres.DrawStim(breakimg);
                        pres.RefreshDisplay;

                    % wait BufferTime 
                        pause(BufferTime)     


                    % EEG processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % get markers and time samples
                        [mark, t_mark] = inlet_marker.pull_chunk(); 
                    % get EEG data and time samples
                        [eeg, eeg_t] = inlet_eeg.pull_chunk();
                    % preprocess trials
                        [Pre_trials, ChExcl] = SD_cleanERPs(eeg, eeg_t, t_mark, TimeInfo);
                    % calculate ERP across channels of interest
                        [ERP_Features, ERPs] = SD_ERPFeatures(Pre_trials, TimeInfo);

                    % checking how many NaN-timeseries across 8 channels X 12 trials 
                         count = sum(ChExcl.Excl(:) == 1)
                         total = size(ChExcl.Excl(:))
                         perc_bad = (count/total(1))*100
                         perc_good = 100 - perc_bad

                    % Channel feedback; included in ERP or not
                         [FigChQual,Incl] = SD_Channel_Feedback_perc(ChExcl, TimeInfo, FigChQual, Ntrials, Blck, N_goodQualERP_aim, perc_good);

                    % Repeat block or continues?
                         prompt = 'Repeat block (r), save block and continue (c) or save block and quit session (q)?';
                         str2 = input(prompt,'s');
                        
                         if  contains(str2,'c')

                     % Report and bookkeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     % plot final ERP with identified features
                                Ntrials = size(Pre_trials,3);
                                [Fig_Nc] = SD_ERPts_Figure(Fig_Nc, ERPs, ERP_Features, Ntrials, Blck, TimeInfo, ' ');

                    % EEGraw data into structure for this block
                                EEGraw = struct();
                                EEGraw.eeg = eeg;
                                EEGraw.eeg_t = eeg_t;
                                EEGraw.mark = mark;
                                EEGraw.t_mark = t_mark;    

                    % append clean trials and latencies
                                if isequal(Blck,1) % for the first block

                    % clean eeg data
                                    ClnTrials = Pre_trials;

                    % Nc
                                    Output.Nc_Lat = ERP_Features.Nc_Lat;
                                    Output.Nc_Ampl = ERP_Features.Nc_Ampl;
                                    Output.Nc_MAmpl = ERP_Features.Nc_MAmpl;

                    % ERP time series
                                    ERPtimeseries{1,1} = ERPs;

                    % Channel and eeg data
                                    Ch_Excl{1,1} = ChExcl;
                                    EEGdata_raw{1,1} = EEGraw;

                                else % concatenate with other blocks
                    % clean eeg data
                                    ClnTrials = cat(4,ClnTrials,Pre_trials); %Ch x Time x Trl x Blck            
                    % Nc
                                    Output.Nc_Lat = cat(1,Output.Nc_Lat,ERP_Features.Nc_Lat);
                                    Output.Nc_Ampl = cat(1,Output.Nc_Ampl,ERP_Features.Nc_Ampl);
                                    Output.Nc_MAmpl = cat(1,Output.Nc_MAmpl,ERP_Features.Nc_MAmpl);
                    % ERP time series    
                                    ERPtimeseries = cat(1,ERPtimeseries,ERPs); %Blck
                    % Channel and eeg data
                                    Ch_Excl = cat(1,Ch_Excl,ChExcl); %Blck
                                    EEGdata_raw = cat(1,EEGdata_raw,EEGraw); %Blck
                                end

                                clear Pre_trials ChExcl 
                                clear ii countPres eeg eeg_t mark t_mark
                                clear ERP_Features

                                NameFig = strcat(Name_pptfolder,'/ERPs_Blck_',num2str(Blck),'.png');
                                saveas(Fig_Nc, NameFig)

                                NameMatBlock = strcat(Name_pptfolder, '/Faces_eeg_output_blck',num2str(Blck),'.mat');    
                                save(NameMatBlock)

                                endblock = sprintf('End of block %2d \n', Blck);
                                disp(endblock);

                                % save variables into .mat file at end of the session
                                Faces = Output;
                                save(NameMat,'ClnTrials','Faces','Ch_Excl', 'EEGdata_raw','TimeInfo','Output')
                                save(strcat(Name_pptfolder,'/EEGdata_raw.mat'), 'EEGdata_raw')


                                Blck = Blck + 1; 
                                i = i + 1;

                         elseif contains(str2,'r')
                             rep_text = sprintf('Repeating block %2d',Blck);
                                    disp(rep_text);
                                    Blck = Blck;
                                    i = i;
                         
                         else 
                     % Report and bookkeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     % plot final ERP with identified features
                                Ntrials = size(Pre_trials,3);
                                [Fig_Nc] = SD_ERPts_Figure(Fig_Nc, ERPs, ERP_Features, Ntrials, Blck, TimeInfo, ' ');

                    % EEGraw data into structure for this block
                                EEGraw = struct();
                                EEGraw.eeg = eeg;
                                EEGraw.eeg_t = eeg_t;
                                EEGraw.mark = mark;
                                EEGraw.t_mark = t_mark;    

                    % append clean trials and latencies
                                if isequal(Blck,1) % for the first block

                    % clean eeg data
                                    ClnTrials = Pre_trials;

                    % Nc
                                    Output.Nc_Lat = ERP_Features.Nc_Lat;
                                    Output.Nc_Ampl = ERP_Features.Nc_Ampl;
                                    Output.Nc_MAmpl = ERP_Features.Nc_MAmpl;

                    % ERP time series
                                    ERPtimeseries{1,1} = ERPs;

                    % Channel and eeg data
                                    Ch_Excl{1,1} = ChExcl;
                                    EEGdata_raw{1,1} = EEGraw;

                                else % concatenate with other blocks
                    % clean eeg data
                                    ClnTrials = cat(4,ClnTrials,Pre_trials); %Ch x Time x Trl x Blck            
                    % Nc
                                    Output.Nc_Lat = cat(1,Output.Nc_Lat,ERP_Features.Nc_Lat);
                                    Output.Nc_Ampl = cat(1,Output.Nc_Ampl,ERP_Features.Nc_Ampl);
                                    Output.Nc_MAmpl = cat(1,Output.Nc_MAmpl,ERP_Features.Nc_MAmpl);
                    % ERP time series    
                                    ERPtimeseries = cat(1,ERPtimeseries,ERPs); %Blck
                    % Channel and eeg data
                                    Ch_Excl = cat(1,Ch_Excl,ChExcl); %Blck
                                    EEGdata_raw = cat(1,EEGdata_raw,EEGraw); %Blck
                                end

                                clear Pre_trials ChExcl 
                                clear ii countPres eeg eeg_t mark t_mark
                                clear ERP_Features

                                NameFig = strcat(Name_pptfolder,'/ERPs_Blck_',num2str(Blck),'.png');
                                saveas(Fig_Nc, NameFig)

                                NameMatBlock = strcat(Name_pptfolder, '/Faces_eeg_output_blck',num2str(Blck),'.mat');    
                                save(NameMatBlock)

                                endblock = sprintf('End of block %2d \n', Blck);
                                disp(endblock);

                                % save variables into .mat file at end of the session
                                Faces = Output;
                                save(NameMat,'ClnTrials','Faces','Ch_Excl', 'EEGdata_raw','TimeInfo','Output')
                                save(strcat(Name_pptfolder,'/EEGdata_raw.mat'), 'EEGdata_raw')

                           return

                         end       


                    end  
                    return
    end  
end
 
 
 
 



% create empty output files    
dlmwrite('/Users/braintools/Desktop/BONDS_SD/Output.txt',[])   
 

close all

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Experiment completed - End of the session')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

