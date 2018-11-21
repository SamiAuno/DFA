%% Detrended Fluctiation Analysis
% Brainstorm process
% Sami Auno 14.11.2018

% note: data on muodossa [channels,time], eli esim [387,72075]

function varargout = process_DFA( varargin )
% PROCESS_DFA: calculates Detrended Fluctuation Analysis for given data
%
% USAGE:      sProcess = process_DFA('GetDescription')
%               sInput = process_DFA('Run',     sProcess, sInput)
%                alpha = process_DFA('Compute', ???)

% ----------------------------------------------------------------------
%
% Author: Sami Auno, University of Helsinki, September 2018
% 
% ----------------------------------------------------------------------

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Detrended Fluctuation Analysis [test]';
    sProcess.FileTag     = 'DFA';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 650;
    sProcess.isSeparator = 1;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'raw', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'raw', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.processDim  = 1;   % Process channel by channel
    sProcess.isSourceAbsolute = 0;
    
    % Definition of the options
    sProcess.options.timewindow.Comment = 'Measurement time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    
    % Split in time blocks
    sProcess.options.split.Comment = 'Split recordings in time blocks (0=disable): ';
    sProcess.options.split.Type    = 'value';
    sProcess.options.split.Value   = {0, 's', []};
    
    sProcess.options.removeSeg.Comment    = 'Remove BAD segments';
    sProcess.options.removeSeg.Type       = 'checkbox';
    sProcess.options.removeSeg.Value      = 1;
    
    sProcess.options.sep.Type    = 'separator';
    sProcess.options.sep.Comment = ' ';
    
    sProcess.options.label2.Comment = '<U><B>DFA analysis paramenters</B></U>';
    sProcess.options.label2.Type    = 'label';
    
    sProcess.options.minTimeWindow.Comment = 'Minimum time window size:';
    sProcess.options.minTimeWindow.Type    = 'value';
    sProcess.options.minTimeWindow.Value   = {0,'s ',3};
    
    sProcess.options.maxTimeWindow.Comment = 'Maximum time window size:';
    sProcess.options.maxTimeWindow.Type    = 'value';
    sProcess.options.maxTimeWindow.Value   = {0,'s ',3};
    
    sProcess.options.numTimeWindows.Comment = 'Number of time windows:';
    sProcess.options.numTimeWindows.Type    = 'value';
    sProcess.options.numTimeWindows.Value   = {20, '',0};
    
    sProcess.options.label1.Comment = '<U><B>Fluctuation calculated by mean or median of the RMS variation?</B></U>';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.avgtype.Comment = {'Mean','Median'};
    sProcess.options.avgtype.Type    = 'radio';
    sProcess.options.avgtype.Value   = 1;
end

% %% ===== GET NUMBER OF TIME SPLITS ===== %%
% function Comment = getNumTimeSplits(sProcess)
% %     Comment = 'Number of blocks:  -';
%     
%     if isfield(sProcess.options, 'split') && isfield(sProcess.options.split, 'Value') && sProcess.options.split.Value{1} > 0 && ~isempty(sProcess.options.timewindow.Value)
%         totTime = sProcess.options.timewindow.Value{2} - sProcess.options.timewindow.Value{1};
% %         nBlocks = totTime/sProcess.options.split{}
%     else
%         Comment = 'Number of blocks:  -';
%     end
% end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    % Get time window
    Comment = [sProcess.Comment, ': [', process_extract_time('GetTimeString', sProcess), ']'];
    % Absolute values 
    if isfield(sProcess.options, 'source_abs') && sProcess.options.source_abs.Value
        Comment = [Comment, ', abs'];
    end
end

%% ===== GET TIME STRING =====
function strTime = GetTimeString(sProcess, sInput)
    % Get time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        time = sProcess.options.timewindow.Value{1};
    elseif (nargin >= 2) && isfield(sInput, 'TimeVector') && ~isempty(sInput.TimeVector)
        time = sInput.TimeVector([1 end]);
    else
        time = [];
    end
    % Print time window
    if ~isempty(time)
        if any(abs(time) > 2)
            if (time(1) == time(2))
                strTime = sprintf('%1.3fs', time(1));
            else
                strTime = sprintf('%1.3fs,%1.3fs', time(1), time(2));
            end
        else
            if (time(1) == time(2))
                strTime = sprintf('%dms', round(time(1)*1000));
            else
                strTime = sprintf('%dms,%dms', round(time(1)*1000), round(time(2)*1000));
            end
        end
    else
        strTime = 'all';
    end
end

%% ===== RUN =====
function sInput = Run(sProcess, sInput)
    % Get time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value{1})
        iTime = panel_time('GetTimeIndices', sInput.TimeVector, sProcess.options.timewindow.Value{1});
    else
        iTime = 1:length(sInput.TimeVector);
    end
    if isempty(iTime)
        bst_report('Error', sProcess, [], 'Invalid time definition.');
        sInput = [];
        return;
    end
    
    if isfield(sProcess.options, 'minTimeWindow') && isfield(sProcess.options.minTimeWindow, 'Value') && iscell(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value{1})
        minTimeWindow = sProcess.options.minTimeWindow.Value{1};
    else
        minTimeWindow = sInput.TimeVector(1);
    end
    if isempty(minTimeWindow)
        bst_report('Error', sProcess, [], 'Invalid minimum time window definition.');
        sInput = [];
        return;
    end
    
    if isfield(sProcess.options, 'maxTimeWindow') && isfield(sProcess.options.maxTimeWindow, 'Value') && iscell(sProcess.options.maxTimeWindow.Value) && ~isempty(sProcess.options.maxTimeWindow.Value) && ~isempty(sProcess.options.maxTimeWindow.Value{1})
        maxTimeWindow = sProcess.options.maxTimeWindow.Value{1}; 
    else
        maxTimeWindow = sInput.TimeVector(end);
    end
    if maxTimeWindow == 0
        maxTimeWindow = sInput.TimeVector(end);
    end
    if isempty(maxTimeWindow)|| maxTimeWindow < minTimeWindow
        bst_report('Error', sProcess, [], 'Invalid maximum time window definition.');
        sInput = [];
        return;
    end
    
    if isfield(sProcess.options, 'numTimeWindows') && isfield(sProcess.options.numTimeWindows, 'Value') && iscell(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value{1})
        numTimeWindows = uint8(sProcess.options.numTimeWindows.Value{1});
    else
        numTimeWindows = uint8(20);
    end
    if isempty(numTimeWindows)
        bst_report('Error', sProcess, [], 'Invalid number of time windows.');
        sInput = [];
        return;
    end
    
    % avgflag for tag
    if sProcess.options.avgtype.Value(1) == 2
        avgtag = 'median';
    elseif sProcess.options.avgtype.Value(1) == 1
        avgtag = 'mean';
    end
    
    data = sInput.A;
    protocolInfo = bst_get('ProtocolInfo');
%     sInput.Fluctuation = [];
    
    % remove segments that contain the word 'remove'
    if (sProcess.options.removeSeg.Value)
        Events = load(fullfile(protocolInfo.STUDIES,sInput.DataFile),'Events');
        Events = Events.Events;

        nDiffEvents = length(Events);

        if(nDiffEvents > 0 && sProcess.options.removeSeg.Value)

            dF = size(sInput.A,2)/(sInput.TimeVector(end) - sInput.TimeVector(1));

            for i=1:nDiffEvents
                if (~isempty(regexp(Events(i).label,'\w*remove', 'once')) || ~isempty(regexp(Events(i).label,'transient', 'once')) )
                    segments = round((Events(i).times(:,:) - sInput.TimeVector(1) ).*dF + 1);
                    if size(segments,1) == 1
                        for j=1:size(segments,2)
                            data(:,segments(1,j)) = nan;
                        end
                    elseif size(segments,1) == 2
                        for j=1:size(segments,2)
                            data(:,segments(1,j):segments(2,j)) = nan;
                        end
                    end
                end
            end

            data = data(:,~isnan(data(1,:)));
            iTime = 1:size(data,2);
        end
    end
    
    frequency = length(sInput.TimeVector)/diff(sInput.TimeVector([1,end]));
    
    % Compute DFA
    [measures,fluctuation] = Compute(data(:, iTime, :),frequency*minTimeWindow,frequency*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));
    sInput.A = measures(:,1);
    
    % save fluctuations
    aux_saveMat(fluctuation, sInput.iBlockRow, fileparts(fullfile(protocolInfo.STUDIES,sInput.DataFile)),'Fluctuations.mat');
    aux_saveMat(measures, sInput.iBlockRow, fileparts(fullfile(protocolInfo.STUDIES,sInput.DataFile)),'Measures.mat');
    % save measures
    
    
    % Copy values to represent the time window
%     sInput.A = [sInput.A, sInput.A];
    % Keep only first and last time values
    sInput.TimeVector = mean([minTimeWindow,sInput.TimeVector(end)]);
%     if (length(iTime) >= 2)
%         sInput.TimeVector = [sInput.TimeVector(iTime(1)), sInput.TimeVector(iTime(end))];
%     % Only one time point: the duplicated time samples must have different time values
%     else
%         if (length(sInput.TimeVector) > 2)
%             sInput.TimeVector = sInput.TimeVector(iTime(1)) + [0, sInput.TimeVector(2)-sInput.TimeVector(1)];
%         else
%             sInput.TimeVector = sInput.TimeVector(iTime(1)) + [0, 1e-6];
%         end
%     end
    % Build file tag
    sInput.CommentTag = [sProcess.FileTag '(' num2str(numTimeWindows) ' win, ' num2str(minTimeWindow,3) '-' num2str(maxTimeWindow,3) ' s, ' avgtag ')'];
    % Do not keep the Std/TFmask fields in the output
    if isfield(sInput, 'Std') && ~isempty(sInput.Std)
        sInput.Std = [];
    end
    if isfield(sInput, 'TFmask') && ~isempty(sInput.TFmask)
        sInput.TFmask = [];
    end
end


%% ===== EXTERNAL CALL =====
function [measures,fluctuation] = Compute(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag)
    
    [measures,fluctuation] = routine_dfa_calc(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag);
    %     [alpha,rsq] = routine_dfa_calc_old(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag);    
    
end














