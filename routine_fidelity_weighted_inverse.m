%% Fidelity weighted inverse solution
% 
% [wIS,msg] = routine_fidelity_weighted_inverse(SourceIdent, FS, IS)
% calculates fidelity weighted inverse solution (wIS).
% SourceIdent is size [n_source,1] vector that contains the parcel index
% for each source.
% FS is the forward solution matrix.
% IS is the inverse solution matrix
%
%
% Based on python implementation: https://github.com/sanrou/fidelityWeighting
% Sami Auno, 2018


function [weighted_inverse_solution,PLVmat,msg] = routine_fidelity_weighted_inverse(source_identities, forward_solution, inverse_solution)
%% First general checking and error throwing.
msg = [];
if isempty(source_identities) 
    msg = 'source_identities vector is empty.';
    return
end
if isempty(forward_solution)
    msg = 'forward_solution matrix is empty.';
    return
end
if isempty(inverse_solution)
    msg = 'inverse_solution matrix is empty.';
    return
end

% Check that the number of vertices matches
if ~(size(source_identities,1) == size(inverse_solution,1) && size(source_identities,1) == size(forward_solution,2))
    msg = 'There is mismatch in number of vertices between source_identities, inverse_solution and forward_solution.';
    return
end

% Check that the number of channels matches.
if ~(size(inverse_solution,2) == size(forward_solution,1))
    msg = 'There is a mismatch in number of channels.';
    return
end

% Get the number of parcels and the number of sources
n_parcels = max(source_identities);
n_sources = size(source_identities,1);

%% Calculate PLV data matrix
% This is saved as a txt file where 
directory = pwd;
PLVdataSaveFileName = 'temp_PLVmat.txt';

% Each time instPLV (instantaneous PLV,) is called, it calculates nDataPoints data points
% To get adequate results, you need minimum total of 100k data points.
% The reason this is done this way is to conserve RAM. If you feel frisky,
% this can be calculated with one step, just increase the number of data
% points in the instPLV function.

% Also, since the PLV data for each data point is saved into text file,
% this can be used for debugging (in which case you need to comment the
% delete() function).

numIterations = 10;
nDataPoints = 10000;
for j=1:numIterations
    disp(['Calculating PLV ', num2str(j), '/', num2str(numIterations)])
    instPLV(source_identities, forward_solution, inverse_solution,nDataPoints,directory,PLVdataSaveFileName)
end

% Read the saved matrix to calculate the nonparametric PLV
instPLVmat = dlmread(fullfile(directory,PLVdataSaveFileName));

% Now that the temp.-file is no longer needed, it can be deleted
delete(fullfile(directory,PLVdataSaveFileName));

% This is the total number of data points
analysisLength = numIterations*nDataPoints;

% And here we calculate the nonparametric PLV for each source
PLVmat = (abs(sum(instPLVmat))/analysisLength)';
%% Calculate weights and create fidelity weighted inverse solution

weights = zeros(n_sources,1);

% The weight is PLV^2. Sign is used as a flip
for i=1:n_sources
    weights(i,1) = PLVmat(i,1)^2 * sign(PLVmat(i,1));
end

% Calculate unnormalized weighted inverse operator
weighted_inverse_solution = sparse(1:n_sources,1:n_sources,weights) * inverse_solution;

% Calculate normalized weights
normalizedWeights = zeros(n_sources,1);
for i=1:n_parcels
    ind_sources = find( source_identities == i );
    
    normalizedWeights(ind_sources,1) = weights(ind_sources,1) * ( norm(inverse_solution(ind_sources,:)) / norm(weighted_inverse_solution(ind_sources,:) ) );
end

% Calculate normalized weighted inverse operator
weighted_inverse_solution = sparse(1:n_sources,1:n_sources,normalizedWeights) * inverse_solution;
weighted_inverse_solution = weighted_inverse_solution * ( norm(inverse_solution) / norm(weighted_inverse_solution) );



end


function instPLV(source_identities, forward_solution, inverse_solution, timePointsFin, directory, saveName)
% This calculates the instantaneous PLV


% Get the number of parcels and the number of sources
n_parcels = max(source_identities);
n_sources = size(source_identities,1);

%% Create random dummy signal for each parcel
% The time series for each parcel is timePointsFin long. Include cutout
% regions with length of timePointsCut to exclude border effects.

timePointsCut = 150;
timePointsGen = timePointsFin + 2*timePointsCut;

parcelTimeSeries = randn(n_parcels,timePointsGen);

% Do convolution with Ricker wavelet
for i=1:n_parcels
    parcelTimeSeries(i,:) = conv(parcelTimeSeries(i,:),mexihat(-5,5,timePointsFin/100),'same');
end

% Do Hilbert transform to get analytic signal.
parcelTimeSeries = hilbert(parcelTimeSeries')';

% Cut borders off
parcelTimeSeries = parcelTimeSeries(:,timePointsCut+1:timePointsGen-timePointsCut);

%% Create dummy source space based on the parcel time series
% This is done by cloning each parcel time series to those sources that are
% within the parcel. The parcel IDs are in the source_identities vector.

sourceTimeSeries = zeros(n_sources,timePointsFin);

for i=1:n_sources
    if source_identities(i,1) == 0 % check for sources that are not part of any parcel
        sourceTimeSeries(i,:) = zeros(1,size(parcelTimeSeries,2));
    else
        sourceTimeSeries(i,:) = parcelTimeSeries(source_identities(i,1),:);
    end
end

%% Create dummy source space based on the source estimation
% Use the forward and inverse models to calculate new source space time
% series. This introduces interparcel signal mixing.

sourceTimeSeries = inverse_solution * ( forward_solution * sourceTimeSeries );

%% Calculate the instantaneous PLV, so PLV as a function of time

PLVt = zeros(n_sources,1);
% PLVt = zeros(n_sources,timePointsFin);

% Change to amplitude 1, maintain angle
% sourceTimeSeries = exp(1j*angle(sourceTimeSeries));
% parcelTimeSeries = exp(1j*angle(parcelTimeSeries));

% relative phase if defined as:
% z1(t) = a(t) + ix(t), z2(t) = b(t) + iy(t)
% R = abs(z), theta = angle(z)
% z = R*exp(i*theta) = R*ẑ, ẑ = unit vector

% PLVt = ẑ1(t) .* conjugate( ẑ2(t) )
% To get PLV independent of time, sum over time, take the magnitude and normalize by number of time points 

for i=1:n_sources
%     PLVt(i,1) = sum( (parcelTimeSeries(source_identities(i,1),:)./(abs(parcelTimeSeries(source_identities(i,1),:))) ) .* conj( (sourceTimeSeries(i,:))./(abs(sourceTimeSeries(i,:))) ) ,2 );
%     PLVt(i,1) = sum( exp(1j*angle(parcelTimeSeries(source_identities(i,1),:))) .* conj( exp(1j*angle(sourceTimeSeries(i,:))) ) ,2 );
    
    if source_identities(i,1) == 0 % check for sources that are not part of any parcel
        PLVt(i,1) = complex(0,0);
    else
        PLVt(i,1) = sum(exp(1j*angle(parcelTimeSeries(source_identities(i,1),:))) .* conj( exp(1j*angle(sourceTimeSeries(i,:))) ),2);
    end
end

% Save to file
dlmwrite(fullfile(directory,saveName),PLVt','-append');

end











