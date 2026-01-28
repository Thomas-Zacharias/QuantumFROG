function [retrievedPulseList, Gt] = svdexFROGMM(fieldSpectrogramListFrogIter, numMode, iterFrog)
% svdexFROG: Extracts the pulse and gate as functions of time from the FROG

% fieldSpectrogramListFrogIter: array of field spectrograms corresponding to different modes
% numMode: number of modes
% iterFrog: current frog iteration
% retrievedPulseList: array of retrieved modes
% Gt: complex gate profile (if using blind frog)

N = size(fieldSpectrogramListFrogIter, 1);

vList = zeros(N, numMode);
uList = zeros(N, numMode);
eList  = zeros(N, numMode); 
projList = zeros(N, numMode); 

modeOccupancyList = zeros(1, numMode);


%% Multi mode Gram Schmidt
for iter = 1:numMode
    N = size(fieldSpectrogramListFrogIter(:, :, iter), 1);
    %Do the exact inverse of the procedure in makeXFROG...
    %Undo the line: EF=circshift(ifft(EF,[],1),ceil(N/2)-1);
    EF = fft(circshift(fieldSpectrogramListFrogIter(:, :, iter), 1-ceil(N/2)), [], 1);
    %Undo the line: EF=fliplr(ifftshift(EF,2));
    EF = fftshift(fliplr(EF),2);
    %Undo the lines: for n=2:N  EF(n,:) = circshift(EF(n,:), [0 1-n]);  end
    for n=2:N
        EF(n,:) = circshift(EF(n,:), [0, (n-1)]);
    end
    % Now EF is the "outer product form", see Kane1999.
    
    try 
        [U, S, V] = svds(EF,1);
    catch
        EF(isnan(EF)) = 0;
        [U, S, V] = svds(EF,1);
    end 

    vList(:, iter) = S.*U(:, 1);    
    if iter == 1
        uList(:, iter) = vList(:, iter);
        eList(:, iter) = uList(:, iter)./(uList(:, iter)'*uList(:, iter));
    else
        % Gram Schmidt
        for iterProj = 1:iter-1
            projList(:, iter) = projList(:, iter)  + (uList(:, iterProj)'*vList(:, iter)).*eList(:, iterProj);
        end
        uList(:, iter) = vList(:, iter) - projList(:, iter);
        eList(:, iter) = uList(:, iter)./(uList(:, iter)'*uList(:, iter));
    end

    Gt = conj(V(:, 1));
    
%     toc 
    %% Adding constraint to reorder from highest to lowest mode occupancy
    if iterFrog > 100
        modeOccupancyList(iter) = norm(uList(:, iter));
    end
end 
retrievedPulseList = uList;

if iterFrog > 100
% Order from highest to lowest mode occupancy
    [~, sortedIndex] = sort(modeOccupancyList, 'descend');
    retrievedPulseList = retrievedPulseList(:, sortedIndex);

end


end