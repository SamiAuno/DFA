%% Detrended Fluctiation Analysis
% Brainstorm process
% Sami Auno

% note: Data matrix should be of size [number of channels,number of timepoints ], eg. [387,72075]

function varargout = process_DFA( varargin )
% PROCESS_DFA: calculates Detrended Fluctuation Analysis for given data
%
% USAGE:      sProcess = process_DFA('GetDescription')
%               sInput = process_DFA('Run',     sProcess, sInput)
%                alpha = process_DFA('Compute')

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
    
%     % Definition of the options
    sProcess.options.timewindow.Comment = 'Measurement time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    
    % Segments (events) that are named either 'transient' (eg. due to
    % filtering) or have the string 'remove' at the end of the name, are
    % not included in the process
    sProcess.options.removeSeg.Comment    = 'Remove transient/bad segments';
    sProcess.options.removeSeg.Type       = 'checkbox';
    sProcess.options.removeSeg.Value      = 1;
    
    sProcess.options.sep.Type    = 'separator';
    sProcess.options.sep.Comment = ' ';
    
    sProcess.options.label2.Comment = '<U><B>DFA analysis paramenters</B></U>';
    sProcess.options.label2.Type    = 'label';
    
    % Minimum time window for DFA
    sProcess.options.minTimeWindow.Comment = 'Minimum DFA time window size:';
    sProcess.options.minTimeWindow.Type    = 'value';
    sProcess.options.minTimeWindow.Value   = {0,'s ',3};
    
    % Maximum time window for DFA
    sProcess.options.maxTimeWindow.Comment = 'Maximum DFA time window size:';
    sProcess.options.maxTimeWindow.Type    = 'value';
    sProcess.options.maxTimeWindow.Value   = {0,'s ',3};
    
    % Number of time windows. The window sizes will be logarithmically
    % spaced between [minTimeWindow,maxTimeWindow]
    sProcess.options.numTimeWindows.Comment = 'Number of DFA time windows:';
    sProcess.options.numTimeWindows.Type    = 'value';
    sProcess.options.numTimeWindows.Value   = {20, '',0};
    
    % The fluctuation may be calculated either with mean or with median.
    % Median is not as susceptible to transient effects
    sProcess.options.label1.Comment = '<U><B>Fluctuation calculated by mean or median of the RMS variation?</B></U>';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.avgtype.Comment = {'Mean','Median'};
    sProcess.options.avgtype.Type    = 'radio';
    sProcess.options.avgtype.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    % Get time window
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== LOAD ALL INFO =====
    % Load the surface filename from results file
    ResultsMat = in_bst_results(sInput.FileName, 0);
    
    % Load variables and chech that there is no conflicts or empty
    % parameters
    
    % Get time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value{1})
        iTime = panel_time('GetTimeIndices', ResultsMat.Time, sProcess.options.timewindow.Value{1});
    else
        iTime = 1:length(ResultsMat.Time);
    end
    
    % Error catching if time empty
    if isempty(iTime)
        bst_report('Error', sProcess, sInput, 'Invalid time definition.');
        return;
    end
    
    % Check that the minimum time window is okay and not conflicting with
    % the maximum time window. Catch errors if needed
    if isfield(sProcess.options, 'minTimeWindow') && isfield(sProcess.options.minTimeWindow, 'Value') && iscell(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value) && ~isempty(sProcess.options.minTimeWindow.Value{1})
        minTimeWindow = sProcess.options.minTimeWindow.Value{1};
    else
        minTimeWindow = 100*(ResultsMat.Time(2)-ResultsMat.Time(1));
    end
    if isempty(minTimeWindow)
        bst_report('Error', sProcess, sInput, 'Invalid minimum time window definition.');
        return;
    end
    
    % Check that the maximum time window is okay and not conflicting with
    % the minimum time window. Catch errors if needed
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
    
    % Check that the number of time windows is okay. Catch errors if
    % needed. If it is empty, then a default 20 windows is used.
    if isfield(sProcess.options, 'numTimeWindows') && isfield(sProcess.options.numTimeWindows, 'Value') && iscell(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value) && ~isempty(sProcess.options.numTimeWindows.Value{1})
        numTimeWindows = uint8(sProcess.options.numTimeWindows.Value{1});
    else
        numTimeWindows = uint8(20);
    end
    if isempty(numTimeWindows) % just in case. You never know.
        bst_report('Error', sProcess, sInput, 'Invalid number of time windows.');
        return;
    end
    
    % Get avgtag. Decided whether to use mean or median to compute the
    % fluctuation
    if sProcess.options.avgtype.Value(1) == 1
        avgtag = 'mean';
    elseif sProcess.options.avgtype.Value(1) == 2
        avgtag = 'median';

    end
    
    % Get protocol info
    protocolInfo = bst_get('ProtocolInfo');
    % Calculate sampling frequency
    samplingFreq = length(ResultsMat.Time)/diff(ResultsMat.Time([1,end]));
    
    % remove events that contain the word 'remove' or 'transient'
    % First all time points are marked as TRUE
    % Then those time points are within the time window of the
    % aforementioned events are marked as FALSE. These are subsequently
    % removed.
    if (sProcess.options.removeSeg.Value)
        
        % Load events list
        Events = load(fullfile(protocolInfo.STUDIES,ResultsMat.DataFile),'Events');
        Events = Events.Events;

        nDiffEvents = length(Events);   % number of events

        if(nDiffEvents > 0)

            iGoodSegments = true(size(ResultsMat.Time));                   % indeces of good segments
            
            for i=1:nDiffEvents
                if (~isempty(regexp(Events(i).label,'\w*remove', 'once')) || ~isempty(regexp(Events(i).label,'transient', 'once')) )
                    segments = round((Events(i).times(:,:) - ResultsMat.Time(1) ).*samplingFreq + 1);
                    if size(segments,1) == 1
                        for j=1:size(segments,2)
                            iGoodSegments(1,segments(1,j)) = false;
                        end
                    elseif size(segments,1) == 2
                        for j=1:size(segments,2)
                            iGoodSegments(1,segments(1,j):segments(2,j)) = false;
                        end
                    end
                end
            end
            iGoodSegments = iGoodSegments(1,iTime);
            iTime = iTime(1,iGoodSegments); % 
        end
    end
    
    
    % ===== COMPUTE DFA =====
    % Do the computation in patches
    n_parcels = size(ResultsMat.ImageGridAmp,1);
    measures = zeros(n_parcels,4);
    fluctuation = zeros(n_parcels,numTimeWindows,3);
    sizeSegs = 20;                          % process 'sizeSegs'-number of channels or parcels at a time
    n_segs = floor((n_parcels-1)/sizeSegs);
    
    for i=1:n_segs
        parc_indices = (1 + sizeSegs*(i-1)):(sizeSegs*i);                   % parcel indices
        [measures(parc_indices,:),fluctuation(parc_indices,:,:)] = Compute(ResultsMat.ImageGridAmp(parc_indices, iTime),samplingFreq*minTimeWindow,samplingFreq*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));
    end
    parc_indices = (1+sizeSegs*n_segs):n_parcels;
    [measures(parc_indices,:),fluctuation(parc_indices,:,:)] = Compute(ResultsMat.ImageGridAmp(parc_indices, iTime),samplingFreq*minTimeWindow,samplingFreq*maxTimeWindow, numTimeWindows, sProcess.options.avgtype.Value(1));
    
    
    
    % ===== SAVE FILE =====
    % Create returned structure 
    NewMat = ResultsMat;
    NewMat.ImageGridAmp  = measures(:,1);   % save alpha as ImageGridAmp
    NewMat.ImagingKernel = [];
    NewMat.Comment       = [NewMat.Comment, ' | DFA (', num2str(numTimeWindows), ' win, ', num2str(minTimeWindow,3), '-', num2str(maxTimeWindow,3), ' s, ', avgtag, ', length ', num2str(round(length(iTime)/samplingFreq)), 's' , ')' ];
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
    
end














