%% Generate theoretical spectrogram/Load measured spectrogram and run multimode retrieval algorithm

%% Initialize
close all; 
clear all;

%% Define units
units = Units();

%% Define signal and pump
signal = Signal();
pump = Pump(); 

%% Define constants
constants = Constants();

%% Theory or experiment
% theory = 1 generates spectrogram from user defined modes and runs retrieval algorithm
% theory = 0 loads measured spectrogram and runs retrieval algorithm
theory = 0;         % Change to 0 to switch to measured spectrogram

%% Load spectrogram
if theory == 1
    % Theoretical Spectrogram: Generates theoretical spectrograms for testing. 
    % Edit TheoreticalSpectram() function to control mode profiles and number of modes
    [t, f, numMode, finalInputSpectrogram, gatePulseEnvelope, intenSpectrogramList, pulseList, squeezingParamList, squeezingLeveldBList] = TheoreticalSpectrogram();        
    [photonNumberFractionalValuesInput, fractionalRetrievedPulseListInput] = PhotonNumberFractionalValues(pulseList, numMode);  % Theoretical Schmidt coefficients 
else
    % ExperimentalSpectrogram(): Loads and preprocess measured spectrogram. 
    % Edit ExperimentalSpectrogram() to change which spectrogram to load
    [t, f, finalInputSpectrogram, gatePulseEnvelope, spectrogramMaxVal] = ExperimentalSpectrogram();
    figure; imagesc(t.axis .* 1e15, f.axis .* 1e-12, finalInputSpectrogram); colorbar;
    xlabel("Time [fs]"); 
    ylabel("Frequency [THz]"); 
    title("Measured spectrogram"); 
    FigureFont();
end

%% Multimode phase retrieval
iterMax = 5000;         % suffuicient for theory grid size: increase iterMax for larger time-grid sizes
% iterMax = 15000;        % for time-grid used in experiments
if theory == 1
    numRetrMode = numMode;   
else
    numRetrMode = 5;    % overestimate and truncate according to retrieved mean photon occupancy for squeezing experiment
%     numRetrMode = 9;    % overestimate and truncate according to retrieved mean photon occupancy for 100THz experiment
end

[retrievedPulseList, retrievedSpectrogramList, retrievedSumLoss, scaledSumRetrievedLoss, iterFrog, retrievedLossList] = svdFROGMM(finalInputSpectrogram, numRetrMode, gatePulseEnvelope, [], [], iterMax, 2, t.step, {'fs','PHz'}, t, f);
[photonNumberFractionalValuesRetr, fractionalRetrievedPulseListRetr] = PhotonNumberFractionalValues(retrievedPulseList, numRetrMode);

figure; plot(photonNumberFractionalValuesRetr); 
title("Schmidt coefficient"); 
xlabel("Modes")
ylabel("Schmidt coefficient")

retrievedSqueezingParamList = zeros(1, numRetrMode);

for iter = 1:numRetrMode
    retrievedSqueezingParamList(iter) = asinh(norm(retrievedPulseList(:, iter)));
end

retrievedSqueezingLeveldBList = 10*log10(exp(2*retrievedSqueezingParamList));

%% Plotting 
if theory == 1
    ReorderAndPlotTheory(t, f, numRetrMode, pulseList, intenSpectrogramList, squeezingParamList, retrievedPulseList, retrievedSpectrogramList, retrievedSumLoss, iterFrog, retrievedLossList, retrievedSqueezingParamList)    
else
    [photonNumberFractionalValues] = NormalizeAndPlotMMExperiment(t, f, numRetrMode, finalInputSpectrogram, retrievedPulseList, retrievedSpectrogramList, retrievedSumLoss, scaledSumRetrievedLoss, iterFrog, retrievedLossList, retrievedSqueezingParamList);    
end

%% Functions
function [t, f, numMode, finalInputSpectrogram, gatePulseEnvelope, intenSpectrogramList, pulseList, squeezingParamList, squeezingLeveldBList] = TheoreticalSpectrogram()
    %% Generate Theoretical spectrogram
    % Edit "numMode" to control number of modes
    % Edit "iterModePulse" to control complex mode profile
    
    units = Units();

    %% Define time axis

    num = 6;   % number of samples = 2^n        % works well with 100fs modes
    timeLength = 1400*units.femto;

    [t, f] = CreateAxis(num, timeLength); 
    
    %% Define FROG Gate pulse
    fwhmPump = sqrt(2)*100*units.femto;    % Defining gate    

    stdev = fwhmPump/(2*sqrt(2*log(2))); 
    gatePulseEnvelope = 1*exp(-(t.axis - t.mid).^2/(2*stdev^2)); 
    gatePulseEnvelope = gatePulseEnvelope./norm(gatePulseEnvelope); % complex field is normalized
    
    %% N HG modes
    numMode = 2;    % Select number of principal modes

    fwhm = sqrt(2)*100*units.femto;

    stdev = fwhm/(2*sqrt(2*log(2))); 
    gdd = (0 * units.femto)^2;
    stdevGdd = sqrt(stdev.^2 + 2i.*gdd);
    stdev = stdevGdd;

    %% Squeezing param to dB level
    squeezingParamList = ones(1, numMode);

    squeezingLeveldBList = 10*log10(exp(2.*squeezingParamList));

    pulseList = zeros(length(t.axis), numMode);
    pulseNormList = zeros(length(t.axis), numMode + 1); % adding vacuum

    figure(30); 

    for iterMode = 1:numMode
        iterModePulse = hermiteH(iterMode-1, (t.axis - t.mid)/stdev).*exp(-(t.axis - t.mid).^2/(2*stdev^2));  
    
        iterOrthoNormModePulse = iterModePulse./norm(iterModePulse);
        pulseNormList(:, iterMode) = iterOrthoNormModePulse;

        pulseList(:, iterMode) = iterOrthoNormModePulse.*sinh(squeezingParamList(iterMode));

        subplot(1, numMode, iterMode); 
        yyaxis left; 
        plot(t.axis./units.femto, abs(iterOrthoNormModePulse)); 
        ylabel("Amplitude [a.u.]");
        yyaxis right; 
        plot(t.axis./units.femto, unwrap(angle(iterOrthoNormModePulse))); 
        ylabel("Phase [rad]");
        title("Mode " + string(iterMode))
    end

    squeezingParamList(numMode+1) = 0;

    %% Make spectrogram
    [intenSpectrogramList, fieldSpectrogramList] = makeXFROGMM(pulseList, gatePulseEnvelope);

    sumIntenSpectrogram = sum(intenSpectrogramList, 3);   
    normVal = 1./max(max(sumIntenSpectrogram));  % testing normalization  
    
    sumIntenSpectrogram = normVal.*sumIntenSpectrogram;
    intenSpectrogramList = normVal.*intenSpectrogramList;
    

    for iterMode = 1:numMode
        figure(57); subplot(numMode+1, 1, iterMode); imagesc(t.axis, f.axis, intenSpectrogramList(:, :, iterMode)); colorbar;
    end
    
    inputSpectrogram = sumIntenSpectrogram;  

    finalInputSpectrogram = zeros(size(inputSpectrogram));
    snrRatio = 1e3;
    finalInputSpectrogram(inputSpectrogram > max(max(inputSpectrogram))/1e3) = inputSpectrogram(inputSpectrogram > max(max(inputSpectrogram))/snrRatio);
    
    figure(57); subplot(numMode+1, 1, numMode+1); imagesc(t.axis, f.axis, finalInputSpectrogram); colorbar;
end

function [t, f, finalSpectrogram, gatePulseEnvelope, spectrogramMaxVal] = ExperimentalSpectrogram()
    %% Generate Experimental spectrogram
    units = Units();

    %% Define Fourier grid    
    num = 9;   % number of samples = 2^n       
    timeLength = 2000*units.femto;

    [t, f] = CreateAxis(num, timeLength); 

    %% Sources
    gate = Pump(); 
    signal = Signal();  

    nonlinearTermCarrierFreq = NonlinearFreq(gate, signal); 

    %% Data - Squeezing chip + SFG FROG
%     fname = 'MeasuredRawData\Data_SqueezingChip\SHGXFROG040125v2.csv';            % amplified vacuum
    fname = 'MeasuredRawData\Data_SqueezingChip\SHGXFROG040125v4ampsq.csv';         % amplified squeezed quadrature
%     fname = 'MeasuredRawData\Data_100THz\XFROG0908v6.csv';                        % 100THz bandwidth
    
    %% Data - OPA + SFG FROG
    if strcmp(fname, 'Data_100THz\XFROG0908v6.csv') == 1
        TiS = 0;
    else
        TiS = 1;
    end

    %% Define FROG Gate pulse
    gatePulseEnvelope = LoadGatePulse(t, TiS);

    %% Import data
    [measPos, measWl, measData1] = ImportData(fname);

    %% Process  measured spectrogram
    %% Crop/Threshold 1
    measData1 = ThresholdBackground(measData1);   % Threshold background noise    
    
    %% Create axes and filter far frequencies
    [measDelay, measFreq, measData1] = CreateAxesFilterFarFrequencies(measPos, measWl, measData1, nonlinearTermCarrierFreq);
 
    %% Threshold 2
    snr = 45;   % for 1mm lihtium iodate
    measData = ThresholdManual(measData1, snr);

    %% Update axes
    [measDelayAxis, measFreqAxis, ac_trace, avespec, carrier_freq, FWHM_ac] = UpdateAxesGenXCorrSpecCorr(measData, measDelay, measFreq);

    %% Normalize and smooth spectrogram
    measData = Smooth(measData);
    
    %% Interpolate spectrogram
    finalSpectrogram = InterpolateSpectrogram(measDelayAxis, measFreqAxis, measData, t, f);
    spectrogramMaxVal= max(max(finalSpectrogram)); 
    finalSpectrogram = finalSpectrogram./spectrogramMaxVal;
end

function measData1 = Smooth(measData)
    filtsig = 0.5; % pixel width of smoothing filter
    measData1 = imgaussfilt(measData,filtsig);
end

function nonlinearTermCarrierFreq = NonlinearFreq(gate, signal) 
    nonlinearTermCarrierFreq = gate.freq + signal.freq; % SFG nonlinearity

    % nonlinearTermCarrierFreq = 4.3103e14; % expected center frequency in Hz, green =
    % 5.747e14, red = 4.3103e14, 1 um = 2.8708e1, split = 5.03e144;
    % exp_carrierfreq2 = 4.0541e14;
    % exp_carrierfreq2 = 4.8387e14; % 620nm
    % nonlinearTermCarrierFreq = 5.747e14;
    % nonlinearTermCarrierFreq = 2.87e14;
    % nonlinearTermCarrierFreq = 5.03e14;
end

function measData1 = ThresholdManual(measData, ratio)
    measData1 = measData;
    measData1(measData1 < max(max(measData1))/ratio) = 0;
    measData1 = measData1 - min(min(measData1(measData1>0)));
    measData1(measData1 < 0) = 0;
end

function interpSpec = InterpolateSpectrogram(measDelayAxis, measFreqAxis, measData, t, f)
    % Rescale to Fourier grid:
    data_timeint = interp1(measDelayAxis, measData, t.axis);   % measdelt_scaled - time we get from stage; measData - corresponding measData; t.axis- our new grid; datatimeint - new measData based on new grid
    interpSpec = interp1(measFreqAxis, data_timeint', f.axis);  % measfreq_scaled - frequency we get from osa; 
    interpSpec(isnan(interpSpec)) = 0;
%     interpSpec = lowpass(interpSpec, 1e-5);
%     interpSpec = (lowpass(interpSpec.', 1e-5)).';
end 

function [measPos, measWl, measData] = ImportData(fname)
    measPos = readmatrix(fname,'Range', '1:1').*1e-3; % position in m
    measWl = readmatrix(fname,'Range', '3:3').*1e-9; % wavelength in m
    
    numtrace = length(measPos);  
    measData = readmatrix(fname,'Range', strcat('4:',num2str(numtrace+3)));
end

function [measDelayAxis, measFreqAxis, ac_trace, avespec, carrier_freq, FWHM_ac] = UpdateAxesGenXCorrSpecCorr(measData, measDelay, measFreq)

    % Find autocorrelation and spectral correlation functions to rescale axes
    ac_trace = mean(measData,2)./max(mean(measData,2));
    ac_trace = ac_trace./max(ac_trace);
    centloc_ac = round(sum((1:length(ac_trace))' .* abs(ac_trace.^4))/sum(abs(ac_trace.^4))); % weighted average to find center
    ac_cent = measDelay(centloc_ac);
    measDelayAxis = -(measDelay - measDelay(centloc_ac)); % center axis at 0 and flip to get the delay axis correct
 
    avespec = mean(measData)./max(mean(measData));
    avespec = avespec./max(avespec);
    centloc_spec = round(sum((1:length(avespec)) .* abs(avespec.^4))/sum(abs(avespec.^4))); % weighted average to find center
    carrier_freq = measFreq(centloc_spec);
    measFreqAxis = (measFreq - measFreq(centloc_spec)); 
    
    % Calculate the FWHM in both domains:
    N_time = length(ac_trace);
    [~,locleft] = min(abs(0.5 - ac_trace(1:centloc_ac)));
    [~,locright] = min(abs(0.5 - ac_trace(centloc_ac:end)));
    locright = locright + (round(centloc_ac)-1);
    FWHM_ac = abs(measDelayAxis(locright) - measDelayAxis(locleft));
end

function [measDelay, measFreq, measData] = CreateAxesFilterFarFrequencies(measpos, measWl, measData, nonlinearTermCarrierFreq)
    constants = Constants();
    BW_coarse = 200*1e12; % Frequency filter bandwidth 
    BW_coarse = 100*1e12; % Frequency filter bandwidth 

    measDelay = measpos*2./constants.c; % convert to temporal delay
    measFreq = constants.c./(measWl); % convert to frequency
    measData(:,abs(measFreq - nonlinearTermCarrierFreq) > BW_coarse) = 0; % filter far frequencies
end

function measData = ThresholdBackground(measData)
    background = measData(1,:);
    measData = measData - background;
    measData(measData < 0) = 0;
end

function gatePulse = LoadGatePulse(t, TiS)
    if TiS == 1
        gatePulse = load("MeasuredRawData\Data_Gates\928nmGateFunction.mat").retrievedPulse;    % different gate
        gateTgrid = load("MeasuredRawData\Data_Gates\928nmGateTimeAxis.mat").t.axis;         
    else
        gatePulse = load("MeasuredRawData\Data_Gates\1umGateFunction.mat").Pt;
        gateTgrid = (load("MeasuredRawData\Data_Gates\1umGateTimeAxis.mat").tgrid)*1e-15; 
    end
    
    gatePulse = interp1(gateTgrid, gatePulse, t.axis);
    gatePulse(isnan(gatePulse)) = 0;    
    gatePulse = gatePulse./norm(gatePulse); 

%     gatePulse = flipud(gatePulse);
end

function [t, f] = CreateAxis(num, timeLength)
    t.N = 2^num;
    t.L = timeLength; 
    t.step = t.L/t.N; 
    t.n = (-t.N/2 + 1:1:t.N/2).';
    t.axis = t.n*t.step; 
    t.length = length(t.axis);
    t.mid = t.axis(round(t.length/2));

    f.N = t.N;
    f.step = 1/t.L;
    f.L = 1/t.step;
    f.n = (-f.N/2:1:f.N/2 - 1).'; 
    f.axis = f.n*f.step;
    f.length = length(f.axis); 
    f.mid = f.axis(round(f.length/2));
end 
     
function [photonNumberFractionalValues, fractionalRetrievedPulseList] = PhotonNumberFractionalValues(fieldList, numMode)
    % photonNumberFractionalValues: Fractional number of photon per mode aka Schmidt Coefficients
    % fractionalRetrievedPulseList: Orthonormal modes weighted by amplitude defined by Schmidt coefficients
    photonNumberTotVal = zeros(1, numMode);
    photonNumberFractionalValues = zeros(1, numMode);
    normRetrievedPulseList = zeros(size(fieldList));
    for iterMode = 1:numMode
        photonNumberTotVal(:, iterMode) = norm(fieldList(:, iterMode)).^2;
    end 
    for iterMode = 1:numMode
        photonNumberFractionalValues(:, iterMode) = norm(fieldList(:, iterMode)).^2./sum(photonNumberTotVal);  
        normRetrievedPulseList(:, iterMode) = fieldList(:, iterMode)./(norm(fieldList(:, iterMode)));
    end 
    fractionalRetrievedPulseList = sqrt(photonNumberFractionalValues).*normRetrievedPulseList;
end

%% Initialziation

function pump = Pump()
%     pump.lambda = 1045*1e-9;

    pump.lambda = 928*1e-9;

    pump.freq = 3*1e8./pump.lambda;   
    pump.period = 1./pump.freq;
    pump.omega = 2*pi.*pump.freq;
end

function signal = Signal()
    signal.lambda = 2090.*1e-9;                 % 2 micron signal wavelength

%     signal.lambda = 1856*1e-9; 

    signal.freq = 3*1e8/signal.lambda;           % 2 micron frequency - ~140 THz
    signal.period = 1./signal.freq;
    signal.omega = 2*pi*signal.freq;         % 2 micron angular frequency ~2*pi*140 radians
end

function [constants] = Constants()
    constants.c = 3*1e8;
end

function [units] = Units()
    units.femto = 1e-15;
    units.pico = 1e-12;
    units.nano = 1e-9;
    units.micro = 1e-6;
    units.milli = 1e-3;
    units.mega = 1e6;
    units.tera = 1e12;
end

%% Plotting functions
function ReorderAndPlotTheory(t, f, numMode, pulseList, intenSpectrogramList, squeezingParamList, retrievedPulseList, retrievedSpectrogramList, retrievedSumLoss, iterFrog, retrievedLossList, retrievedSqueezingParamList)
    %% N mode plotting - Theoretical - known modes
    reorderList = 1:numMode;
    
    reorderRetrievedPulseList = retrievedPulseList;
    reorderRetrievedSpectrogramList = retrievedSpectrogramList;
    reorderRetrievedLossList = zeros(size(retrievedLossList));
    reorderRetrievedSqueezingParamList = zeros(size(retrievedSqueezingParamList));
    
    for iter = 1:numMode
        reorderRetrievedPulseList(:, iter) = retrievedPulseList(:, reorderList(iter));    
        reorderRetrievedSpectrogramList(:, :, iter) = retrievedSpectrogramList(:, :, reorderList(iter));    
        reorderRetrievedLossList(iter) = retrievedLossList(reorderList(iter));
        reorderRetrievedSqueezingParamList(iter) = retrievedSqueezingParamList(reorderList(iter));
    end
    
    PlotMMTheory(t, f, pulseList, intenSpectrogramList, reorderRetrievedSpectrogramList, reorderRetrievedPulseList, reorderRetrievedLossList, retrievedSumLoss, iterFrog, squeezingParamList, reorderRetrievedSqueezingParamList); 
end 

function FigureFont()
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15)
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');

    legendObj = findobj(gcf, 'Type', 'Legend');
    for iterLegend = 1:length(legendObj)
        legendObj(iterLegend, 1).FontSize = 9;
    end

%     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2);
%     pbaspect([1, 1, 1]);
end

function [photonNumberFractionalValues] = NormalizeAndPlotMMExperiment(t, f, numMode, finalInputSpectrogram, retrievedPulseList, retrievedSpectrogramList, retrievedSumLoss, scaledSumLoss, iterFrog, retrievedLossList, retrievedSqueezingParamList)    

    %% Normalization section - normalize such that total photon number = 1
    photonNumberTotVal = zeros(1, numMode);
    photonNumberFractionalValues = zeros(1, numMode);
    normRetrievedPulseList = zeros(size(retrievedPulseList));
    for iterMode = 1:numMode
        photonNumberTotVal(:, iterMode) = norm(retrievedPulseList(:, iterMode)).^2;
    end 
    for iterMode = 1:numMode
        photonNumberFractionalValues(:, iterMode) = norm(retrievedPulseList(:, iterMode)).^2./sum(photonNumberTotVal);  
        normRetrievedPulseList(:, iterMode) = retrievedPulseList(:, iterMode)./(norm(retrievedPulseList(:, iterMode)));
    end 
    fractionalRetrievedPulseList = sqrt(photonNumberFractionalValues).*normRetrievedPulseList;

    %% Spectrogram panel
    figure(); hold on; 
        
    xlimVec = [-250, 250]; % 1mm BBO
    xlimVec = [-450, 450]; % 50um BBO try 1
    
    ylimVec = [-30, 30];  % 1mm BBO
    ylimVec = [-55, 55];  % 50um BBO try 1

    iterCol = 1;
    for iterPlot = 1:(2*(numMode + 1))
        subplot(2, (numMode+1), iterPlot)
        if mod(iterPlot, (numMode + 1)) == 1    % save first column in each row for sum spectrogram
            continue
        else
            if floor(iterPlot/(numMode + 2)) == 0   % row 1
                imagesc(t.axis .* 1e15, f.axis./1e12, retrievedSpectrogramList(:, :, iterPlot-1)); colorbar;
%                 title("Ret mode " + string(iterPlot-1)); 
                title("Retrieved Mode " + string(iterPlot-1) + " Spectrogram");
                xlabel('Time [fs]'); 
                ylabel('Frequency [THz]');    
%                 xlim([-250, 250]);
%                 ylim([-30, 30]);

                xlim(xlimVec);
                ylim(ylimVec);

%                 axis square;
            else                                    % row 2
                yyaxis left; 
                plot(t.axis.*1e15, abs(fractionalRetrievedPulseList(:, iterCol)));
                ylabel("Amplitude [a.u]");
                xlabel("Time [fs]");
                yyaxis right; 
                plot(t.axis.*1e15, angle(fractionalRetrievedPulseList(:, iterCol)));
                xlim(xlimVec);

                ylabel("Phase [rad]"); 
%                 title("N = " + string(retrievedSqueezingParamList(iterCol)));
                title("N = " + string(photonNumberFractionalValues(:, iterCol)));
                legend(["Amplitude", "Phase"])
                title("Retrieved Mode " + string(iterCol));
                iterCol = iterCol + 1;
%                 axis square;
            end
        end 
    end

    subplot(2, (numMode+1), 1); 
    imagesc(t.axis .* 1e15, f.axis./1e12, finalInputSpectrogram); colorbar;
    title("Measured spectrogram + iter = " + string(iterFrog)); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');
    xlim(xlimVec);
    ylim(ylimVec);
%     axis square;
    subplot(2, (numMode+1), (numMode+2)); 
    imagesc(t.axis .* 1e15, f.axis./1e12, sum(retrievedSpectrogramList, 3)); colorbar;
    title("Retrieved/Scaled loss = " + string(retrievedSumLoss) + "/" + string(scaledSumLoss)); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');
    xlim(xlimVec);
    ylim(ylimVec);
%     axis square;

    
    annotation(gcf,'textbox', [0.075, 0.75, 0.020, 0.028], 'String',{'(a)'}, 'Edgecolor', 'none');
    annotation(gcf,'textbox', [0.075, 0.25, 0.020, 0.028], 'String',{'(b)'}, 'Edgecolor', 'none');
%     ax = gcf;
%     exportgraphics(ax,'Experiment.pdf')

    FigureFont();
end

function PlotMMTheory(t, f, pulseList, inputSpectrogramList, retrievedSpectrogramList, retrievedPulseList, retrievedLossList, retrievedSumLoss, iterFrog, squeezingParamList, retrievedSqueezingParamList) 

    %% Spectrogram panel
    figure(); hold on; 
        
    numMode = size(retrievedSpectrogramList, 3);

    tempList = zeros(size(retrievedSpectrogramList));    
    tempList(:, :, 1:size(inputSpectrogramList, 3)) = inputSpectrogramList;
    inputSpectrogramList = tempList;

    tempPulseList = zeros(size(retrievedPulseList));
    tempPulseList(:, 1:size(pulseList, 2)) = pulseList;
    pulseList = tempPulseList;

    tempSqList = zeros(size(retrievedSqueezingParamList));
    tempSqList(:, 1:size(squeezingParamList, 2)) = squeezingParamList; 
    squeezingParamList = tempSqList;

    iterCol = 1;
    for iterPlot = 1:(2*(numMode + 1))
        subplot(3, (numMode+1), iterPlot)
        if mod(iterPlot, (numMode + 1)) == 0
            continue
        else
            if floor(iterPlot/(numMode + 1)) == 0
                imagesc(t.axis .* 1e15, f.axis./1e12, inputSpectrogramList(:, :, iterPlot)); colorbar;
%                 title("In " + string(iterPlot)); 
                ylim([-30 30]);
                title("Input Mode " + string(iterPlot)); 
                xlabel('Time [fs]'); 
                ylabel('Frequency [THz]');        
            else
                imagesc(t.axis .* 1e15, f.axis./1e12, retrievedSpectrogramList(:, :, iterCol)); colorbar;
%                 title("Ret " + string(iterCol) + " " + string(retrievedLossList(iterCol))); 
                ylim([-30 30]);
                title("Retrieved Mode " + string(iterCol)); 
                xlabel('Time [fs]'); 
                ylabel('Frequency [THz]');
                iterCol = iterCol + 1;
            end
        end 
    end

    subplot(3, (numMode+1), numMode+1); 
    imagesc(t.axis .* 1e15, f.axis./1e12, sum(inputSpectrogramList, 3)); colorbar;
%     title("Original FROG iter = " + string(iterFrog)); 
    ylim([-30 30]);
    title("Multimode Spectrogram"); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');

    subplot(3, (numMode+1), 2*(numMode+1)); 
    imagesc(t.axis .* 1e15, f.axis./1e12, sum(retrievedSpectrogramList, 3)); colorbar;
    ylim([-30 30]);
%     title("Sum Loss = " + string(retrievedSumLoss)); 
    title("Retrieved Spectrogram"); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');


    %% Plot panel 
    for iterPlot = 1:numMode
        subplot(3, numMode+1, (2*(numMode+1) + iterPlot))
        yyaxis left; 
        hold on; 
        plot(t.axis.*1e15, abs(retrievedPulseList(:, iterPlot)./norm(retrievedPulseList(:, iterPlot))));
        plot(t.axis.*1e15, abs(pulseList(:, iterPlot))./norm(pulseList(:, iterPlot)));
        ylabel("Amplitude [a.u]");
        xlabel("Time [fs]");
        yyaxis right; 
        plot(t.axis.*1e15, angle(retrievedPulseList(:, iterPlot)));
        plot(t.axis.*1e15, angle(pulseList(:, iterPlot)));
        ylabel("Phase [rad]"); 
%         title("Og/Ret = " + string(retrievedSqueezingParamList(iterPlot) + "/" + string(squeezingParamList(iterPlot))));        
        title("Mode " + string(iterPlot));        
%         legend(["Retrieved amplitude", "Original amplitude", "Retrieved phase", "Original phase"]);
        if t.L > 1e-12
            xlim([-500, 500])
        end
%         xlim([-500, 500])
    end
    
    annotation(gcf,'textbox', [0.075, 0.80, 0.020, 0.028], 'String',{'(a)'}, 'Edgecolor', 'none');
    annotation(gcf,'textbox', [0.075, 0.48, 0.020, 0.028], 'String',{'(b)'}, 'Edgecolor', 'none');
    annotation(gcf,'textbox', [0.075, 0.16, 0.020, 0.028], 'String',{'(c)'}, 'Edgecolor', 'none');

    FigureFont();
end
