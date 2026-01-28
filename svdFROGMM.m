function [retrievedPulseMatFrogIter, retrievedSpectrogramList, sumRetrievedLoss, scaledSumRetrievedLoss, iterFrog, retrievedLossListFrogIter] = svdFROG(inputSpectrogram, numMode, gatePulseEnvelope, seed, GTol, iterMAX, mov, dtperpx, units, t, f)  
    % svdFROGMM: Phase retrieval - apply data constraint and nonlinearity constraint
    
    %   Returned variables: 
    %   retrievedPulseMatFrogIter:          TxN array       Retrieved N principal complex modes of size T
    %   sumRetrievedLoss:                   scalar          RMS loss between input multimode spectrogram and recovered multimode spectrogram 
    %   scaledSumRetrievedLoss:             scalar          RMS loss after scaling recovered spectrogram to best match input spectrogram, see DeLong1996
    %   iterFrog:                           scalar          number of FROG iterations
    %   retrievedLossListFrogIter           scalar          RMS loss per mode between data constraint and retrieved spectrograms 
    
    
    %   inputSpectrogram:                   TxT array       Input 2D spectrogram 
    %   numMode                             scalar          Number of modes to be retrieved
    %   gatePulseEnvelope                   Tx1 array       Complex gate envelope
    %   seed                                Tx1 array       Seed - initial guess. If empty, will start with gaussian noise
    %   GTol                                scalar          Error threshold
    %   iterMAX                             scalar          Maximum number of iterations allowed
    %   mov	                                scalar          0 (default): No updates while solving. 1: Output movie. 2: Print text.
    %	dtperpx                             scalar          Difference in time-delay (dt) between consecutive pixels (default = 1). This automatically fixes frequency units too. Doesn't affect algorithm, only plots and such.
    %   units                               Cell-array      units{1} is units of dtperpx, units{2} is units{1}^-1, the units of frequency. Default: {'a.u.','1/a.u.'} (a.u. means arbitrary units). A common input here would be {'fs','PHz'}
    %   t                                   Tx1 array       time axis 
    %   f                                   Tx1 array       frequency axis
    
    
    rmsdiff = @(F1, F2) sqrt(mean(mean((F1-F2).^2))); %RMS difference in the entries of two real matrices/vectors.
    normalizemax1 = @(M) M/max(max(M)); %normalize a matrix or vector for its maximum to be 1. Must have real nonnegative entries.
    calcalpha = @(Fm,Fr) sum(sum(Fm.*Fr))/sum(sum(Fr.^2)); %calculates alpha, the positive number that minimizes rmsdiff(Fm,alpha*Fr). See DeLong1996
    
    %   Get trace dimensions
    N = size(inputSpectrogram, 1);
    
    % Default behavior for various inputs...
    if (~exist('iterMAX', 'var')||isempty(iterMAX))
        iterMAX = inf;
    end
    if (~exist('GTol', 'var')||isempty(GTol))
        GTol = 0;
    end
    if (~exist('mov', 'var')||isempty(mov))
        mov = 0;
    end
    if (~exist('dtperpx','var')||isempty(dtperpx))
	    dtperpx = 1;
    end
    
    dvperpx = 1 / (N*dtperpx); %frequency interval per pixel
    
    tpxls=(-dtperpx*(N-1)/2:dtperpx:dtperpx*(N-1)/2)'; %x-axis labels for plots
    vpxls=(-dvperpx*(N-1)/2:dvperpx:dvperpx*(N-1)/2)'; %x-axis labels for plots
    
    %maybe you only want to display part of the plot range, to zoom in on the
    %interesting stuff. If so, edit the following lines...
    tplotrange = [min(tpxls) max(tpxls)];
    vplotrange = [min(vpxls) max(vpxls)];
    
    
    %   ------------------------------------------------------------
    %   S T A R T   O F   A L G O R I T H M
    %   ------------------------------------------------------------    
    
    %% Multimode phase retrieval algorithm using Separable State Generalized Projections Algorithm (MSGPA + Gram Schmidt)
    
    
    retrievedPulseMatFrogIter = zeros(size(gatePulseEnvelope, 1), numMode);
    
    if (~exist('seed', 'var')||isempty(seed))
        for iterTemp = 1:numMode
    %             retrievedPulseMatFrogIter(:, iterTemp) = exp(-2*log(2)*(((0:N-1)'-N/2)/(N/20)).^2).*exp(0.1*2*pi*1i*rand(N,1));
            retrievedPulseMatFrogIter(:, iterTemp) = rand(N,1).*exp(0.1*2*pi*1i*rand(N,1));
        end
    end 
    
    if (~exist('gatePulseEnvelope', 'var')||isempty(gatePulseEnvelope))
        Gt = exp(-2*log(2)*(((0:N-1)'-N/2)/(N/10)).^2).*exp(0.1*2*pi*1i*rand(N,1));
    else
        Gt = gatePulseEnvelope;
    end
        
    % Generate FROG trace
    iterFrog = 0;
    
    %EFr is reconstructed FROG trace complex amplitudes ( Fr=|EFr|.^2 )
    [intenSpectrogramListFrogIter, fieldSpectrogramListFrogIter] = makeXFROGMM(retrievedPulseMatFrogIter,Gt);
    
    retrievedLossListFrogIter = zeros(numMode, 1);
    alphaValList = zeros(numMode, 1);
    dataConstraintListFrogIter = zeros(size(gatePulseEnvelope, 1), size(gatePulseEnvelope, 1), numMode);
    
    sumIntenSpectrogramFrogIter = sum(intenSpectrogramListFrogIter, 3); 
    
    for iterTemp = 1:numMode
        dataConstraintListFrogIter(:, :, iterTemp) = sqrt(intenSpectrogramListFrogIter(:, :, iterTemp).*inputSpectrogram./sumIntenSpectrogramFrogIter);
    %         alphaValList(iterTemp) = calcalpha(dataConstraintListFrogIter(:, :, iterTemp).^2, intenSpectrogramListFrogIter(:, :, iterTemp));
    %         intenSpectrogramListFrogIter(:, :, iterTemp) = intenSpectrogramListFrogIter(:, :, iterTemp).*alphaValList(iterTemp);
        retrievedLossListFrogIter(iterTemp) = rmsdiff(dataConstraintListFrogIter(:, :, iterTemp), intenSpectrogramListFrogIter(:, :, iterTemp));
    end
    
    sumRetrievedLoss = rmsdiff(inputSpectrogram, sumIntenSpectrogramFrogIter);
    
    maxNormInputSpectrogram = inputSpectrogram./max(max(inputSpectrogram));
    scaledSpectrogram = sumIntenSpectrogramFrogIter * calcalpha(maxNormInputSpectrogram,sumIntenSpectrogramFrogIter); %scale Fr to best match Fm, see DeLong1996
    scaledSumRetrievedLoss = rmsdiff(maxNormInputSpectrogram,scaledSpectrogram);
    
    
    %% Debug
    %     PlotPulsesGate(t, retrievedPulseList, Gt, iterFrog);
    %     PlotSpectrograms(t, f, intenSpectrogramListFrogIter, dataConstraintListFrogIter, sumIntenSpectrogramFrogIter, inputSpectrogram, retrievedLossListFrogIter, sumRetrievedLoss, iterFrog)
    
    %   ------------------------------------------------------------
    %   F R O G   I T E R A T I O N   A L G O R I T H M
    %   ------------------------------------------------------------
    plotPer = 100;
    
    while (iterFrog<iterMAX) && scaledSumRetrievedLoss > GTol
        iterFrog = iterFrog+1;                  %   keep count of no. of iterations
        if mov==2 && mod(iterFrog, plotPer) == 0
    %             disp(['Iteration number: ' num2str(iterFrog) '  L1: ' num2str(retrievedLossListFrogIter(1)) ' L2 ' num2str(retrievedLossListFrogIter(2)) ' L ' num2str(sumRetrievedLoss)]);
            disp(['Iteration number: ' num2str(iterFrog) ' Abs loss ' num2str(sumRetrievedLoss) ' Scaled loss ' num2str(scaledSumRetrievedLoss)]);
        end
        
        %% Data constraint
    %         tic
        for iterTemp = 1:numMode
            intenSpectrogramListFrogIter(intenSpectrogramListFrogIter(:, :, iterTemp) == 0) = NaN;
    
            dataConstraintListFrogIter(:, :, iterTemp) = sqrt(intenSpectrogramListFrogIter(:, :, iterTemp).*inputSpectrogram./sumIntenSpectrogramFrogIter);
            
            fieldSpectrogramListFrogIter(:, :, iterTemp) = fieldSpectrogramListFrogIter(:, :, iterTemp).*dataConstraintListFrogIter(:, :, iterTemp)./abs(fieldSpectrogramListFrogIter(:, :, iterTemp));
         
            fieldSpectrogramListFrogIter(isnan(fieldSpectrogramListFrogIter(:, :, iterTemp))) = 0;        %   Remove divide by zeros
        end
    %         toc 
    
        %% Nonlinearity constraint
    %         tic
        
        [retrievedPulseMatFrogIter,~] = svdexFROGMM(fieldSpectrogramListFrogIter, numMode, iterFrog);       %   Extract pulse field from FROG complex amplitude
        
    %         toc 
        
        %% Preparing data constraint
    %         tic
    
        [intenSpectrogramListFrogIter, fieldSpectrogramListFrogIter] = makeXFROGMM(retrievedPulseMatFrogIter,Gt);    %   Make a FROG trace from new fields
        
    %         toc 
    
        %% Storing useful values for debugging
    %         tic 
    
        for iterTemp = 1:numMode
            dataConstraintListFrogIter(isnan(dataConstraintListFrogIter(:, :, iterTemp))) = 0;
            intenSpectrogramListFrogIter(isnan(intenSpectrogramListFrogIter(:, :, iterTemp))) = 0;
    
    %             alphaValList(iterTemp) = calcalpha(dataConstraintListFrogIter(:, :, iterTemp).^2, intenSpectrogramListFrogIter(:, :, iterTemp));
    %             intenSpectrogramListFrogIter(:, :, iterTemp) = alphaValList(iterTemp).*intenSpectrogramListFrogIter(:, :, iterTemp);
            retrievedLossListFrogIter(iterTemp) = rmsdiff(dataConstraintListFrogIter(:, :, iterTemp), intenSpectrogramListFrogIter(:, :, iterTemp));
        end
    
        sumIntenSpectrogramFrogIter = sum(intenSpectrogramListFrogIter, 3);
        sumRetrievedLoss = rmsdiff(inputSpectrogram, sumIntenSpectrogramFrogIter);
    
    %         PlotPulsesGate(t, retrievedPulseMatFrogIter, Gt, iterFrog);
    %         PlotSpectrograms(t, f, intenSpectrogramListFrogIter, dataConstraintListFrogIter, sumIntenSpectrogramFrogIter, inputSpectrogram, retrievedLossListFrogIter, sumRetrievedLoss, iterFrog)
       
        maxNormInputSpectrogram = inputSpectrogram./max(max(inputSpectrogram));
        scaledSpectrogram = sumIntenSpectrogramFrogIter * calcalpha(maxNormInputSpectrogram,sumIntenSpectrogramFrogIter); %scale Fr to best match Fm, see DeLong1996
        
    %         if rmsdiff(maxNormInputSpectrogram,scaledSpectrogram) > scaledSumRetrievedLoss      % if loss in current iteration > loss in previous iteration
    %             iterMAX = 1;
    %         end
        
        scaledSumRetrievedLoss = rmsdiff(maxNormInputSpectrogram,scaledSpectrogram);
        
        if mod(iterFrog, plotPer) == 0 || iterFrog == 1
            PlotPulsesGate(t, retrievedPulseMatFrogIter, Gt, iterFrog);
    %             exportgraphics(gcf, 'MMEvolPulse.gif', 'Append', true)
            PlotSpectrograms(t, f, intenSpectrogramListFrogIter, dataConstraintListFrogIter, sumIntenSpectrogramFrogIter, inputSpectrogram, retrievedLossListFrogIter, sumRetrievedLoss, iterFrog)
    %             exportgraphics(gcf, 'MMEvolSpectrogram.gif', 'Append', true)
            
        end
    
    %         toc
    
    end
    %pause(.001)
    %     PlotPulsesGate(t, retrievedPulseList, Gt, iterFrog);
    %     PlotSpectrograms(t, f, intenSpectrogramListFrogIter, dataConstraintListFrogIter, sumIntenSpectrogramFrogIter, inputSpectrogram, retrievedLossListFrogIter, sumRetrievedLoss, iterFrog)
    
    retrievedSpectrogramList = intenSpectrogramListFrogIter;
    for iterMode = 1:numMode
        retrievedPulseMatFrogIter(:, iterMode) = retrievedPulseMatFrogIter(:, iterMode).*exp(-1i*angle(retrievedPulseMatFrogIter(t.N/2 + 3, iterMode)));
    end
    %   ------------------------------------------------------------
    %   E N D   O F   A L G O R I T H M
    %   ------------------------------------------------------------    

end

function PlotPulsesGate(t, retrievedPulseList, Gt, iterFrog)
    pulseListLength = size(retrievedPulseList, 2);

    figure(70);  
    for iterPlot = 1:pulseListLength
        subplot(1, pulseListLength + 1, iterPlot); 
        yyaxis left; 
        plot(t.axis, abs(retrievedPulseList(:, iterPlot))); 
        ylabel("Amplitude [a.u]");
        xlabel("Time");
        yyaxis right; 
        plot(t.axis, angle(retrievedPulseList(:, iterPlot))); 
        ylabel("Phase [rad]"); 
        title("Pulse 1 " + string(iterFrog))
    %     xlim([-40*1e-15 40*1e-15])
    %     xlim([-100*1e-15 100*1e-15])
    end
    
    subplot(1, pulseListLength + 1, pulseListLength + 1); 
    yyaxis left; 
    plot(t.axis, abs(Gt)); 
    ylabel("Amplitude [a.u]");
    xlabel("Time");
    yyaxis right; 
    plot(t.axis, angle(Gt)); 
    ylabel("Phase [rad]"); 
    title("Gate")
end

function PlotSpectrograms(t, f, intenSpectrogramListFrogIter, dataConstraintList, sumIntenSpectrogramFrogIter, inputSpectrogram, retrievedLossList, sumRetrievedLoss, iterFrog)
    figure(71);  

    numMode = size(intenSpectrogramListFrogIter, 3);
    iterRow = 1;
    for iterPlot = 1:(2*numMode)
        subplot(numMode + 1, 2, iterPlot);
        if mod(iterPlot, 2) ~= 0
            imagesc(t.axis .* 1e15, f.axis./1e12, intenSpectrogramListFrogIter(:, :, iterRow)); colorbar;
            title("Spectrogram " + string(iterRow) + " "+ string(retrievedLossList(iterRow))); 
            xlabel('Time [fs]'); 
            ylabel('Frequency [THz]');        
        else
            imagesc(t.axis .* 1e15, f.axis./1e12, dataConstraintList(:, :, iterRow).^2); colorbar;
            title("Data constraint " + string(iterRow)); 
            xlabel('Time [fs]'); 
            ylabel('Frequency [THz]');
            iterRow = iterRow + 1;
        end
    end  
    subplot(numMode + 1, 2, 2*(numMode + 1) -1); 
    imagesc(t.axis .* 1e15, f.axis./1e12, sumIntenSpectrogramFrogIter); colorbar;
    title("Spectrogram sum " + string(sumRetrievedLoss)); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');
    subplot(numMode + 1, 2, 2*(numMode + 1)); 
    imagesc(t.axis .* 1e15, f.axis./1e12, inputSpectrogram); colorbar;
    title("Spectrogram sum " + string(sumRetrievedLoss)); 
    xlabel('Time [fs]'); 
    ylabel('Frequency [THz]');
end