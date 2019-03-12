% this function computes the data fluctuation
function [F,F_var] = routine_fluctuation(DATA, n, N_data, N_chan, avgFlag)

    div = floor(2*(N_data - n)/n) + 1;

    % with 2-dimensional matrix
    % divide data into segments of length of n with 50% overlap
    DATA_segments = zeros(N_chan*div,n);
    detrended_DATA_segments = zeros(N_chan*div,n);
    detrended_DATA_std = zeros(N_chan,div);

    % reshape data
    for di=1:div
        DATA_segments((di-1)*N_chan+1:di*N_chan,:) = DATA(:,((di-1)*n/2+1):((di-1)*n/2+n));
    end

    % Detrend segments
    for d=1:N_chan*div
        detrended_DATA_segments(d,:) = detrend(DATA_segments(d,:));
    end

    % calculate the standard deviation
    detrended_DATA_segments_std = std(detrended_DATA_segments,0,2);

    % rereshape data
    for di=1:div
        detrended_DATA_std(:,di) = detrended_DATA_segments_std(((di-1)*N_chan+1):(di*N_chan), 1);
    end
    
    % three sigma exclusion threshold
    % take the standard deviation of detrended_DATA_std over all channels
    % and exclude all values that are over 3 sigmas away from the median
%     F_std = std(detrended_DATA_std(:));
%     Fmedian = median(detrended_DATA_std(:));
%     exclusion_threshold = [Fmedian - 3*F_std, Fmedian + 3*F_std];
    
    F_var = var(detrended_DATA_std,0,2);
    
    % exclude outliers and calculate fluctuation channel by channel
%     counts = zeros(N_chan,1);
    F = zeros(N_chan,1);
    for ch=1:N_chan
        y = detrended_DATA_std(ch,:);
%         y = y(y>exclusion_threshold(1) & y<exclusion_threshold(2));
        
%         counts(ch,1) = length(y)/div;
        if avgFlag == 2
            F(ch,1) = median(y,2);
        elseif avgFlag == 1 || isempty(avgFlag)
            F(ch,1) = mean(y,2);
        end
    end
    
end