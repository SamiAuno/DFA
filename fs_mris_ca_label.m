function message = fs_mris_ca_label(fsdir)

% init return value;
message = [];

% Freesurfer requires linux or mac
if(~isunix)
    message = 'Requires a UNIX system.';
    return
end

currDir = fileparts(mfilename('fullpath'));

% Check if the custom_atlases folder exists
if ~exist(fullfile(currDir,'custom_atlases'),'dir')
    message = 'Directory "custom_atlases" does not exist in the /process folder.';
    return
end

% Check that the required Freesurfer files exist
if ( exist(fullfile(fsdir,'label'),'dir') && exist(fullfile(fsdir,'surf'),'dir') )
    if ~exist(fullfile(fsdir,'surf','lh.sphere.reg'),'file')
        message = 'File "lh.sphere.reg" does not exist in the Freesurfer /surf folder.';
        return
    elseif ~exist(fullfile(fsdir,'surf','rh.sphere.reg'),'file')
        message = 'File "rh.sphere.reg" does not exist in the Freesurfer /surf folder.';
        return
    end
else
    message = 'Directories "/label" or "/surf" do not exist in the Freesurfer subject folder.';
    return
end

% Get atlas directory
atlasesDir = fullfile(currDir,'custom_atlases');
% Get subjects directory and the subject name
[subjsDir,subjName,~] = fileparts(fsdir);
% Get subject canonical surface files
lh_canonSurfFile = fullfile(fsdir,'surf','lh.sphere.reg');
rh_canonSurfFile = fullfile(fsdir,'surf','rh.sphere.reg');

% Open folder selection dialog
atlasDir = java_getfile( 'open', ...
    'Select new atlas(es) to import...', ...     % Window title
    atlasesDir, ...           % custom_atlases directory
    'multiple', 'dirs', ...                  % Selection mode
    {{'.folder'}, 'FreeSurfer folder', 'atlasDir'}, 0);
% If no folder was selected: exit
if isempty(atlasDir)
    return
end

Nnew_atlases = size(atlasDir,1);

% Loop over each new atlas
for ik = 1:Nnew_atlases
    % extract the atlas name
    [~,atlas,~] = fileparts(atlasDir{ik});
    
    % get classifier paths and define output filename
    switch atlas
        case 'hcp-mmp'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.hcp-mmp_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.hcp-mmp_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.hcpmmp.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.hcpmmp.annot');
            
        case 'schaefer100-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer100-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer100-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer100.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer100.annot');
            
        case 'schaefer200-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer200-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer200-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer200.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer200.annot'); 
        case 'schaefer400-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer400-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer400-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer400.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer400.annot');
        case 'schaefer600-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer600-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer600-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer600.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer600.annot');
        case 'schaefer800-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer800-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer800-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer800.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer800.annot');
        case 'schaefer1000-yeo17'
            lh_atlas_classifier = fullfile(atlasDir{ik},'lh.schaefer1000-yeo17_6p0.gcs');
            rh_atlas_classifier = fullfile(atlasDir{ik},'rh.schaefer1000-yeo17_6p0.gcs');
            
            lh_outputfile = fullfile(fsdir,'label','lh.schaefer1000.annot');
            rh_outputfile = fullfile(fsdir,'label','rh.schaefer1000.annot');
        otherwise
            message = 'Unknown atlas.';
            return
    end
    
    % Parse freesurfer commands
    lh_command = ['mris_ca_label -sdir ', subjsDir, ' ',subjName, ' lh ', lh_canonSurfFile, ' ', lh_atlas_classifier, ' ', lh_outputfile];
    rh_command = ['mris_ca_label -sdir ', subjsDir, ' ' ,subjName, ' rh ', rh_canonSurfFile, ' ', rh_atlas_classifier, ' ', rh_outputfile];
    
    
    % Execute freesurfer commands
    % first left hemisphere
    % status = 0, everything is ok
    
    fprintf('Start relabeling surfaces... \n\n');
    status_lh = system(lh_command,'-echo');
    status_rh = system(rh_command,'-echo');
    if ~(status_lh == 0 && status_rh == 0)
        message = 'Errors in freesurfer execution. Exit function.';
        return
    end
    fprintf('Relabeling done \n\n');
end




end
%left hemi	: mris_ca_label -sdir ./subjects/ 3306 lh ./subjects/3306/surf/lh.sphere.reg ./atlas_data/hcp-mmp/lh.hcp-mmp_6p0.gcs ./subjects/3306/label/lh.hcpmmp.annot
%right hemi : mris_ca_label -sdir ./subjects/ 3306 rh ./subjects/3306/surf/rh.sphere.reg ./atlas_data/hcp-mmp/rh.hcp-mmp_6p0.gcs ./subjects/3306/label/rh.hcpmmp.annot
