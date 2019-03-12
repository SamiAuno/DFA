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
% Author: Sami Auno, University of Helsinki, 2018
% 
% ----------------------------------------------------------------------

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Detrended Fluctuation Analysis [test]';
    sProcess.FileTag     = 'DFA';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 650;
    sProcess.isSeparator = 1;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
%     sProcess.processDim  = 1;   % Process channel by channel
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
    
    sProcess.options.minTimeWindow.Comment = 'Minimum DFA time window size:';
    sProcess.options.minTimeWindow.Type    = 'value';
    sProcess.options.minTimeWindow.Value   = {0,'s ',3};
    
    sProcess.options.maxTimeWindow.Comment = 'Maximum DFA time window size:';
    sProcess.options.maxTimeWindow.Type    = 'value';
    sProcess.options.maxTimeWindow.Value   = {0,'s ',3};
    
    sProcess.options.numTimeWindows.Comment = 'Number of DFA time windows:';
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
    Comment = sProcess.Comment;
end

% %% ===== GET TIME STRING =====
% function strTime = GetTimeString(sProcess, sInput)
%     % Get time window
%     if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
%         time = sProcess.options.timewindow.Value{1};
%     elseif (nargin >= 2) && isfield(sInput, 'TimeVector') && ~isempty(sInput.TimeVector)
%         time = sInput.TimeVector([1 end]);
%     else
%         time = [];
%     end
%     % Print time window
%     if ~isempty(time)
%         if any(abs(time) > 2)
%             if (time(1) == time(2))
%                 strTime = sprintf('%1.3fs', time(1));
%             else
%                 strTime = sprintf('%1.3fs,%1.3fs', time(1), time(2));
%             end
%         else
%             if (time(1) == time(2))
%                 strTime = sprintf('%dms', round(time(1)*1000));
%             else
%                 strTime = sprintf('%dms,%dms', round(time(1)*1000), round(time(2)*1000));
%             end
%         end
%     else
%         strTime = 'all';
%     end
% end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== LOAD ALL INFO =====
    % Load the surface filename from results file
    ResultsMat = in_bst_results(sInput.FileName, 0);
    % Get time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value{1})
        iTime = panel_time('GetTimeIndices', ResultsMat.Time, sProcess.options.timewindow.Value{1});
    else
        iTime = 1:length(ResultsMat.Time);
    end
    if isempty(iTime)
        bst_report('Error', sProcess, sInput, 'Invalid time definition.');
        return;
    end
    
    if isfield(sProcess.options, 'minTimeWindow') && isfield(sProcess.options.minTimeWindow, 'Value') && iscell(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value{1})
        minTimeWindow = sProcess.options.minTimeWindow.Value{1};
    else
        minTimeWindow = 100*(ResultsMat.Time(2)-ResultsMat.Time(1));
    end
    if isempty(minTimeWindow)
        bst_report('Error', sProcess, sInput, 'Invalid minimum time window definition.');
        return;
    end
    
    if isfield(sProcess.options, 'maxTimeWindow') && isfield(sProcess.options.maxTimeWindow, 'Value') && iscell(sProcess.options.maxTimeWindow.Value) && ~isempty(sProcess.options.maxTimeWindow.Value) && ~isempty(sProcess.options.maxTimeWindow.Value{1})
        maxTimeWindow = sProcess.options.maxTimeWindow.Value{1}; 
    else
        maxTimeWindow = (ResultsMat.Time(iTime(end))-ResultsMat.Time(iTime(1)))/10; % 10 % of the recording length
    end
    if maxTimeWindow == 0
        maxTimeWindow = (ResultsMat.Time(iTime(end))-ResultsMat.Time(iTime(1)))/10; % 10 % of the recording length
    end
    if isempty(maxTimeWindow)|| maxTimeWindow < minTimeWindow
        bst_report('Error', sProcess, sInput, 'Invalid maximum time window definition.');
        return;
    end
    
    if isfield(sProcess.options, 'numTimeWindows') && isfield(sProcess.options.numTimeWindows, 'Value') && iscell(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value{1})
        numTimeWindows = uint8(sProcess.options.numTimeWindows.Value{1});
    else
        numTimeWindows = uint8(20);
    end
    if isempty(numTimeWindows)
        bst_report('Error', sProcess, sInput, 'Invalid number of time windows.');
        return;
    end
    
    % avgflag for tag
    if sProcess.options.avgtype.Value(1) == 2
        avgtag = 'median';
    elseif sProcess.options.avgtype.Value(1) == 1
        avgtag = 'mean';
    end
    
%     data = sInput.A;
    protocolInfo = bst_get('ProtocolInfo');
    freq = length(ResultsMat.Time)/diff(ResultsMat.Time([1,end]));
%     sInput.Fluctuation = [];
    
    % remove segments that contain the word 'remove'
    if (sProcess.options.removeSeg.Value)
        Events = load(fullfile(protocolInfo.STUDIES,ResultsMat.DataFile),'Events');
        Events = Events.Events;

        nDiffEvents = length(Events);

        if(nDiffEvents > 0)

%             dF = size(sInput.A,2)/(sInput.TimeVector(end) - sInput.TimeVector(1));
            iGoodSegments = logical(ones(size(ResultsMat.Time)));                   % indeces of good segments
            
            for i=1:nDiffEvents
                if (~isempty(regexp(Events(i).label,'\w*remove', 'once')) || ~isempty(regexp(Events(i).label,'transient', 'once')) )
                    segments = round((Events(i).times(:,:) - ResultsMat.Time(1) ).*freq + 1);
                    if size(segments,1) == 1
                        for j=1:size(segments,2)
                            iGoodSegments(1,segments(1,j)) = false;
%                             data(:,segments(1,j)) = nan;
                        end
                    elseif size(segments,1) == 2
                        for j=1:size(segments,2)
                            iGoodSegments(1,segments(1,j):segments(2,j)) = false;
%                             data(:,segments(1,j):segments(2,j)) = nan;
                        end
                    end
                end
            end
            iGoodSegments = iGoodSegments(1,iTime);
%             data = data(:,~isnan(data(1,:)));
            iTime = iTime(1,iGoodSegments);
        end
    end
    
%     frequency = length(ResultsMat.Time)/diff(ResultsMat.Time([1,end]));
   

%     % Check the length of the data and whether it is over 10% of the
%     % maximum time window. If not, change the maximum timewindow to be 10%
%     % of the lenght of the resulting measurement
%     
%     if 0.1*length(ResultsMat.ImageGridAmp(1, iTime)) < freq*maxTimeWindow
%         fprintf('process_DFA> maxTimeWindow over 10 %% of data length. \n');
%         fprintf('process_DFA> Changing maxTimeWindow to be 10 %% of data length \n\n');
%         maxTimeWindow = (0.1*length(ResultsMat.ImageGridAmp(1, iTime)))/freq;
%     end
    
    % ===== COMPUTE DFA =====
    % Do the computation in patches
    n_parcels = size(ResultsMat.ImageGridAmp,1);
    measures = zeros(n_parcels,4);
    fluctuation = zeros(n_parcels,numTimeWindows,3);
    sizeSegs = 20;                          % process 'sizeSegs'-number of channels or parcels at a time
    n_segs = floor((n_parcels-1)/sizeSegs);
    
    for i=1:n_segs
        parc_indices = (1 + sizeSegs*(i-1)):(sizeSegs*i);                   % parcel indices
        [measures(parc_indices,:),fluctuation(parc_indices,:,:)] = Compute(ResultsMat.ImageGridAmp(parc_indices, iTime),freq*minTimeWindow,freq*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));
    end
    parc_indices = (1+sizeSegs*n_segs):n_parcels;
    [measures(parc_indices,:),fluctuation(parc_indices,:,:)] = Compute(ResultsMat.ImageGridAmp(parc_indices, iTime),freq*minTimeWindow,freq*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));
    
    
    
    
    %     [measures,fluctuation] = Compute(ResultsMat.ImageGridAmp(:, iTime),freq*minTimeWindow,freq*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));

    
%     sInput.A = measures(:,1);
    
    % save fluctuations
%     aux_saveMat(fluctuation, sInput.iBlockRow, fileparts(fullfile(protocolInfo.STUDIES,sInput.DataFile)),'Fluctuations.mat');
%     aux_saveMat(measures, sInput.iBlockRow, fileparts(fullfile(protocolInfo.STUDIES,sInput.DataFile)),'Measures.mat');
    % save measures
    
    
%     % Copy values to represent the time window
% %     sInput.A = [sInput.A, sInput.A];
%     % Keep only first and last time values
%     sInput.TimeVector = mean([minTimeWindow,sInput.TimeVector(end)]);
% %     if (length(iTime) >= 2)
% %         sInput.TimeVector = [sInput.TimeVector(iTime(1)), sInput.TimeVector(iTime(end))];
% %     % Only one time point: the duplicated time samples must have different time values
% %     else
% %         if (length(sInput.TimeVector) > 2)
% %             sInput.TimeVector = sInput.TimeVector(iTime(1)) + [0, sInput.TimeVector(2)-sInput.TimeVector(1)];
% %         else
% %             sInput.TimeVector = sInput.TimeVector(iTime(1)) + [0, 1e-6];
% %         end
% %     end
%     % Build file tag
%     sInput.CommentTag = [sProcess.FileTag '(' num2str(numTimeWindows) ' win, ' num2str(minTimeWindow,3) '-' num2str(maxTimeWindow,3) ' s, ' avgtag ')'];
%     % Do not keep the Std/TFmask fields in the output
%     if isfield(sInput, 'Std') && ~isempty(sInput.Std)
%         sInput.Std = [];
%     end
%     if isfield(sInput, 'TFmask') && ~isempty(sInput.TFmask)
%         sInput.TFmask = [];
%     end
    % ===== SAVE FILE =====
    % Create returned structure 
    NewMat = ResultsMat;
    NewMat.ImageGridAmp  = measures(:,1);   % save alpha as ImageGridAmp
    NewMat.ImagingKernel = [];
    % NewMat.Comment       = [NewMat.Comment, ' | atlas' num2str(length(sScouts))];
    NewMat.Comment       = [NewMat.Comment, ' | DFA (', num2str(numTimeWindows), ' win, ', num2str(minTimeWindow,3), '-', num2str(maxTimeWindow,3), ' s, ', avgtag, ', length ', num2str(round(length(iTime)/freq)), 's' , ')' ];
    NewMat.Time          = [];
    NewMat.Measures      = measures;
    NewMat.Fluctuation   = fluctuation;
    NewMat.AnalysisLength = length(iTime);
    % Add history entry
    NewMat = bst_history('add', NewMat, 'dfa', ['DFA: ', num2str(numTimeWindows), ' win, ', num2str(minTimeWindow,3), '-', num2str(maxTimeWindow,3), ' s, ', avgtag]);
    % Get output study
    sStudy = bst_get('Study', sInput.iStudy);
    % File tag
    fileTag = 'results_dfa';
    % Output filename
    OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), fileTag);
    % Save on disk
    bst_save(OutputFiles{1}, NewMat, 'v6');
    % Register in database
    db_add_data(sInput.iStudy, OutputFiles{1}, NewMat);
end


%% ===== EXTERNAL CALL =====
function [measures,fluctuation] = Compute(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag)
    
    [measures,fluctuation] = routine_dfa_calc(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag);
    %     [alpha,rsq] = routine_dfa_calc_old(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag);    
    
end














