function [OutputVal] = SD_ERPoutputVal(ERP_Features, Outputmeasure)

% The function [OutputVal] = SD_ERPoutputVal(ERP_Features, Outputmeasure)
% selects the output measure of interest for Bayesian Optimisation.

% input:
% - ERP_Features; ERP measures
% - Outputmeasure; measure of interest

% Edits: RH 12-09-19, adapted by ET 08-11-19, AG 02-12-19, AG 4-01-21

%% Select output measure and create report string

switch Outputmeasure 
    case 'Nc_Ampl' 
            formatSpec = 'Nc Amplitude: %.2f uV';
            OutputVal = abs(ERP_Features.Nc_Ampl);  
    case 'Nc_Lat' 
            formatSpec = 'Nc Latency: %.2f ms';
            OutputVal = ERP_Features.Nc_Lat;
    case 'AUC_Nc_window' 
            formatSpec = 'AUC for Nc window: %.2f ';
            OutputVal = ERP_Features.AUC_Nc_window;   
    case 'Nc_MAmpl' 
            formatSpec = 'MAmplitude Nc window: %.2f uV';
            OutputVal = ERP_Features.Nc_MAmpl*(-1);  
            
end

    % save output value
    dlmwrite('/Users/braintools/Desktop/BONDS_SD/Output.txt',OutputVal) 
    
 
    
end
