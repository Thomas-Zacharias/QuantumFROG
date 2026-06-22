%% Load amplified vacuum and amplified squeezed vacuum data post phase retrieval alogrithm and recovers quantum pulse
%% Generates figures used in Quantum FROG paper

close all; 
clear all;

%% Units
units = Units();

%% Amplified vacuum 
% Load recovered data
load("RecoveredModesData\ampVac.mat");
PlotSpectrogramsAmpPhaseDataOrthonormAndGainAndSpectrograms(ampVac, ampVac.orthonormcomplexModes, ampVac.gainMat, ampVac.complexModes, ampVac.gainMat, " "); 

%% Amplified squeeezed vacuum 
% Load recovered data
load("RecoveredModesData\ampSqVac.mat");

% Normalizing factor for measured spectrogram integration time
sqVacRelNormFactor = 0.6137;            
ampSqVac.measuredSpectrogram = sqVacRelNormFactor * ampSqVac.measuredSpectrogram;

PlotSpectrogramsAmpPhaseDataOrthonormAndGainAndSpectrograms(ampSqVac, ampSqVac.orthonormcomplexModes, ampSqVac.gainMat, ampSqVac.complexModes, ampVac.gainMat, " "); 

%% Plot temmporal correlation matrix
PlotTempCorrMat(ampVac.t.axis, ampVac.tempCorrMat, ampSqVac.tempCorrMat);

%% Invert action of measurement amplifier
CharacterizeQuantumPulse(ampVac, ampSqVac);

%% Functions 
function PlotTempCorrMat(tAxis, ampVacTempCorrMat, ampSqVacTempCorrMat)
    %% Plot temporal correlation matrix
    units = Units(); 
    figure; 
    subplot(1, 4, 1); 
    imagesc(tAxis./units.femto, tAxis./units.femto, abs(ampVacTempCorrMat)./max(max(abs(ampVacTempCorrMat))));colorbar; 
    ylabel("Time [fs]"); 
    xlabel("Time [fs]"); 
    ylim([-150 150]); 
    xlim([-150, 150]);
    title("Amp vac temp corr mat")
    axis square; 

    subplot(1, 4, 2); 
    imagesc(tAxis./units.femto, tAxis./units.femto, angle(ampVacTempCorrMat));colorbar; 
    ylabel("Time [fs]"); 
    xlabel("Time [fs]"); 
    ylim([-150 150]); 
    xlim([-150, 150]); 
    title("Phase vac temp corr mat")
    axis square; 

    subplot(1, 4, 3); 
    imagesc(tAxis./units.femto, tAxis./units.femto, abs(ampSqVacTempCorrMat)./max(max(abs(ampVacTempCorrMat))));colorbar; 
    ylabel("Time [fs]"); 
    xlabel("Time [fs]"); 
    ylim([-150 150]); 
    xlim([-150, 150]);
    title("Amp sq vac temp corr mat")
    axis square; 
%     clim([0, max(max(abs(ampVacTempCorrMat)))])
    clim([0, 1])

    subplot(1, 4, 4); 
    imagesc(tAxis./units.femto, tAxis./units.femto, angle(ampSqVacTempCorrMat));colorbar; 
    ylabel("Time [fs]"); 
    xlabel("Time [fs]"); 
    ylim([-150 150]); 
    xlim([-150, 150]); 
    title("Phase Sq vac temp corr mat")
    axis square; 
    WhiteColorMap();

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

function PlotSpectrogramsAmpPhaseDataOrthonormAndGainAndSpectrograms(ampQuad, orthonormcomplexModes, gainMat, complexModes, ampVacGainMat, sqStr)
    %% Plot measured and recovered spectrogram, and recovered complex modes and associated spectrograms

    units = Units();

    TiS = 1;
    gatePulseEnvelope = LoadGatePulse(ampQuad.t, TiS);
    [intenSpectrogramList, fieldSpectrogramList] = makeXFROGMM(complexModes, gatePulseEnvelope);

    figure; 
    imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, squeeze(sum(intenSpectrogramList, 3)));
    xlim([-400, 400]); 
    ylim([-40, 40]); 
    xlabel("Time [fs]"); 
    ylabel("Frequency [THz]")
    axis square; 

    colorbar; clim([0, 1])
    WhiteColorMap();

    for iter = 1:size(orthonormcomplexModes, 2)
        subplot(2, size(orthonormcomplexModes, 2) + 1, iter+1); 
        yyaxis left; 
        plot(ampQuad.t.axis./units.femto, abs(orthonormcomplexModes(:, iter))); 
        xlabel("Time [fs]"); 
        ylabel("Amplitude [a.u.]"); 
        xlim([-150 150]); 

        yyaxis right; 
        plot(ampQuad.t.axis./units.femto, angle(orthonormcomplexModes(:, iter))); 
        xlabel("Time [fs]"); 
        ylabel("Phase [rad.]"); 
        xlim([-200 200]); 

        if strcmp(sqStr, "Squeezing") == 1
            title("Mode " + string(iter) + " Squeezing " + string(10*log10(gainMat(iter, iter).^2)) + " dB")
        else 
            title("Mode " + string(iter) + " gain " + string(gainMat(iter, iter).^2 ./ sum(diag(ampVacGainMat).^2)))
        end

        axis square; 

        subplot(2, size(orthonormcomplexModes, 2) + 1, size(orthonormcomplexModes, 2) + iter+2)
        imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, squeeze(intenSpectrogramList(:, :, iter)));
        xlim([-400, 400]); 
        ylim([-40, 40]); 
        xlabel("Time [fs]"); 
        ylabel("Frequency [THz]")
        axis square; 

        colorbar; 
    end

    subplot(2, size(orthonormcomplexModes, 2) + 1, 1); 
    imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, ampQuad.measuredSpectrogram);
    xlim([-400, 400]); 
    ylim([-40, 40]); 
    xlabel("Time [fs]"); 
    ylabel("Frequency [THz]")
    axis square; 
    title("Measured spectrogram")
    colorbar;

    subplot(2, size(orthonormcomplexModes, 2) + 1, size(orthonormcomplexModes, 2) + 2); 
    imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, squeeze(sum(intenSpectrogramList, 3)));
    xlim([-400, 400]); 
    ylim([-40, 40]); 
    xlabel("Time [fs]"); 
    ylabel("Frequency [THz]")
    axis square; 
    title("Recovered spectrogram")

    WhiteColorMap();

%     subplot(1, size(orthonormcomplexModes, 2) + 1, iter + 1);
%     imagesc(gainMat); colorbar;   axis square; title("Gain")
end

function PlotAmpPhaseDataOrthonormAndGainAndSpectrograms(ampQuad, orthonormcomplexModes, gainMat, complexModes, ampVacGainMat, sqStr)
    %% Plot recovered complex modes and associated spectrograms
    units = Units(); 

    TiS = 1;
    gatePulseEnvelope = LoadGatePulse(ampQuad.t, TiS);
    [intenSpectrogramList, fieldSpectrogramList] = makeXFROGMM(complexModes, gatePulseEnvelope);

    figure; 
    imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, squeeze(sum(intenSpectrogramList, 3)));
    xlim([-400, 400]); 
    ylim([-40, 40]); 
    xlabel("Time [fs]"); 
    ylabel("Frequency [THz]")
    axis square; 

    colorbar; clim([0, 1])
    WhiteColorMap();

    for iter = 1:size(orthonormcomplexModes, 2)
        subplot(2, size(orthonormcomplexModes, 2), iter); 
        yyaxis left; 
        plot(ampQuad.t.axis./units.femto, abs(orthonormcomplexModes(:, iter))); 
        xlabel("Time [fs]"); 
        ylabel("Amplitude [a.u.]"); 
        xlim([-150 150]); 

        yyaxis right; 
        plot(ampQuad.t.axis./units.femto, angle(orthonormcomplexModes(:, iter))); 
        xlabel("Time [fs]"); 
        ylabel("Phase [rad.]"); 
        xlim([-200 200]); 

        if strcmp(sqStr, "Squeezing") == 1
            title("Mode " + string(iter) + " Squeezing " + string(10*log10(gainMat(iter, iter).^2)) + " dB")
        else 
            title("Mode " + string(iter) + " gain " + string(gainMat(iter, iter).^2 ./ sum(diag(ampVacGainMat).^2)))
        end

        axis square; 

        subplot(2, size(orthonormcomplexModes, 2), iter + size(orthonormcomplexModes, 2))
        imagesc(ampQuad.t.axis./units.femto, ampQuad.f.axis./units.tera, squeeze(intenSpectrogramList(:, :, iter)));
        xlim([-400, 400]); 
        ylim([-40, 40]); 
        xlabel("Time [fs]"); 
        ylabel("Frequency [THz]")
        axis square; 

        colorbar; 
    end

    WhiteColorMap();

end

function CharacterizeQuantumPulse(ampVac, ampSqVac)
    units = Units(); 

    %% Find temporal correlation matrix of amplified vacuum 
    ampVac.tempCorrMat = ampVac.complexModes * ctranspose(ampVac.complexModes);

    %% Find temporal correlation matrix of ampSqVac
    ampSqVac.tempCorrMat = ampSqVac.complexModes * ctranspose(ampSqVac.complexModes);    % \Psi\dagger_{TxM}\Psi_{MxT}
    
    %% Inverting measurement OPA 
    ampSqVac.inputModeCorrMat = abs(inv(ampVac.gainMat)  * ctranspose(ampVac.orthonormcomplexModes) * ampSqVac.tempCorrMat * ampVac.orthonormcomplexModes * inv(ampVac.gainMat));   % Inverting measurement amplifier and applying physical constraint corresponding to high gain amplifier approximation
    
    figure;
    subplot(1, 2, 1)
    imagesc(abs(ampSqVac.inputModeCorrMat)); colorbar; axis square; 
    title("Amp: In sq quad corr mat"); xlabel("Modes"); ylabel("Modes"); 
    subplot(1, 2, 2); 
    imagesc(ampSqVac.t.axis./units.femto, ampSqVac.t.axis./units.femto, angle(ampSqVac.inputModeCorrMat)); axis square; 
    title("Phase: In sq quad corr mat"); colorbar; xlabel("Modes"); ylabel("Modes");
    
    figure; 
    imagesc(ampSqVac.inputModeCorrMat); colorbar; axis square; 
    clim([0, 1])
    WhiteColorMap();

    %% Diagonalize
    [ampSqVac.eigenVecMat, ampSqVac.eigenValMat] = eig(ampSqVac.inputModeCorrMat);
    figure;
    subplot(1, 3, 1)
    imagesc(abs(ampSqVac.eigenVecMat)); colorbar; axis square; 
    title("Amp: eigenmat")
    subplot(1, 3, 2); 
    imagesc(angle(ampSqVac.eigenVecMat)); axis square; 
    title("Phase: eigenmat"); colorbar; 
    subplot(1, 3, 3); 
    imagesc(abs(ampSqVac.eigenValMat));axis square; 
    title("Eigenval"); colorbar; 

    %% Temporal mode transformation
    ampVac.inputOrthonormComplexModes = conj(ampVac.orthonormcomplexModes);

    ampSqVac.inputTempCorrMat =  ampVac.inputOrthonormComplexModes * ampSqVac.inputModeCorrMat * ctranspose(ampVac.inputOrthonormComplexModes);
    [eigMatVec, eigMatVal] = eig(ampSqVac.inputTempCorrMat);
        
    ampSqVac.inputTempModeMat = eigMatVec(:, 1:size(ampSqVac.inputModeCorrMat, 1));
    ampSqVac.inputEigVal = eigMatVal(1:size(ampSqVac.inputModeCorrMat, 1), 1:size(ampSqVac.inputModeCorrMat, 1));
    ampSqVac.inputGainMat = sqrt(ampSqVac.inputEigVal);
    ampSqVac.inputComplexModes = ampSqVac.inputTempModeMat * ampSqVac.inputGainMat;

    PlotAmpPhaseDataOrthonormAndGainAndSpectrograms(ampVac, ampSqVac.inputTempModeMat, ampSqVac.inputGainMat, ampSqVac.inputComplexModes, ampVac.gainMat, "Squeezing"); 
    
end

function WhiteColorMap()
    
    n = 256;                          % total number of colors
    fade_fraction = 0.2;              % fraction that should fade from white
    nfade = round(n * fade_fraction);
    
    % Original jet colormap
    jetmap = jet(n - nfade);
    
    % White first jet color gradient
    fade = [linspace(1, jetmap(1,1), nfade)', ...
            linspace(1, jetmap(1,2), nfade)', ...
            linspace(1, jetmap(1,3), nfade)'];
    
    % Concatenate smooth fade with the rest of jet
    cmap = [fade; jetmap];
    
    colormap(cmap);
    colorbar;

%     cmap = jet(n);
    blend_strength = 0.3;   % 0 = original jet, 1 = all white
    cmap = (1-blend_strength)*cmap + blend_strength*[1 1 1];
    
    colormap(cmap);
    colorbar;
    
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
    