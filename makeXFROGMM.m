function [intenSpectrogramList, fieldSpectrogramList] = makeXFROGMM(pulseList, gatePulseEnvelope)
% makeXFROGMM: Generates SFG-XFROG multimode spectrogram (sum of single mode spectrograms)

% pulseList: Vector containing complex modes
% gatePulseEnvelope: complex Gate envelope used in XFROG measurement
% intenSpectrogramList: array of intensity spectrograms correspoding to different modes
% fieldSpectrogramList: array of field spectrograms corresponding to different modes

N = size(pulseList, 1);

fieldSpectrogramList = zeros(N, N, size(pulseList, 2)); 
intenSpectrogramList = zeros(N, N, size(pulseList, 2)); 

for iter = 1:size(pulseList, 2)    
    fieldSpectrogramList(:, :, iter) = pulseList(:, iter)*(gatePulseEnvelope.');

    %   Row rotation...eqns (10)-->(11) of Kane1999
    for n=2:N
	    fieldSpectrogramList(n,:, iter) = circshift(fieldSpectrogramList(n,:, iter), [0 1-n 0]);
    end
    
    % EF is eqn (11) of Kane1999. From left column to right column, it's
    % tau=0,-1,-2...3,2,1
    
    %permute the columns to the right order, tau=...,-1,0,1,...
    fieldSpectrogramList(:, :, iter) = fliplr(ifftshift(fieldSpectrogramList(:, :, iter),2));
    %FFT each column and put 0 frequency in the correct place:
    fieldSpectrogramList(:, :, iter) = circshift(ifft(fieldSpectrogramList(:, :, iter),[],1),ceil(N/2)-1);
        
    % Generate FROG trace (= |field|^2)

%     fieldSpectrogramList(:, :, iter) = fieldSpectrogramList(:, :, iter)./max(max(fieldSpectrogramList(:, :, iter)));    % normalization 
    
    intenSpectrogramList(:, :, iter) = abs(fieldSpectrogramList(:, :, iter)).^2;
end