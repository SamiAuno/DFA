%% Calculate DFA exponent
% measures = [alpha,rsq,b_index,conf]
% Sami Auno, 14.11.2018


function [measures,fluctuation] = routine_dfa_calc(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag)

% if minSampleWindow < 10
%     minSampleWindow = 10;
% end

[N_chan,N_data] = size(DATA);
fluctuation = nan(N_chan,numTimeWindows,3);

% calcuale Hilbert transform -> envelope
fprintf('routine_dfa_calc> Calculating signal envelope... \n\n');
DATA_uenv = zeros(N_chan,N_data);
for i=1:N_chan
    [DATA_uenv(i,:),~] = envelope(DATA(i,:));
end
clear DATA;
fprintf('routine_dfa_calc> Signal envelope calculated. \n\n');

% number of time windows in DFA
res = numTimeWindows;                   

% initializing and defining variables

n_vec = 2*floor(logspace(log10(minSampleWindow),log10(maxSampleWindow),res)/2);   % only even numbers, spread logarithmically

% Check that the longest time window (max(n_vec)) is not over 10% of the
% total data length. In case it is, exclude those time windows that are
% over 10% from the data set
if n_vec(end) > N_data*0.1
    fprintf('routine_dfa_calc> maxTimeWindow over 10 %% of data length. \n');
    
    n_vec = n_vec(n_vec <=N_data*0.1);
    res = uint8(length(n_vec));
    
    fprintf('routine_dfa_calc> Only consides the first %d/%d windows, which are under 10 %% of data length. \n', res,numTimeWindows);
    
end


F_n = zeros(N_chan,res); % fluctuation as a function of n_vec
F_var = zeros(N_chan,res);
alpha = zeros(N_chan,1);
rsq = zeros(N_chan,1);
b_index = zeros(N_chan,1);
conf = zeros(N_chan,1);

% compute cumulative sum.
DATA_cs = zeros(N_chan,N_data);
for i=1:N_chan
    DATA_cs(i,:) = cumsum(DATA_uenv(i,:),2);
end
clear DATA_uenv;

fprintf('routine_dfa_calc> Estimating fluctuation function... \n\n');

% compute fluctuation for each n_vec
for j=1:res

    [F_n(:,j),F_var(:,j)] = routine_fluctuation(DATA_cs, n_vec(j), N_data, N_chan, avgFlag);

end
clear DATA_cs;
fprintf('routine_dfa_calc> Fluctuation functions calculated. \n\n');

fprintf('routine_dfa_calc> Estimating DFA exponents... \n\n');

% linear regression analysis, compute alpha
log10F_n = log10(F_n);
n_vec = log10(n_vec);

fluctuation(:,1:res,1) = log10F_n;
fluctuation(:,1:res,2) = repmat(n_vec,[N_chan,1]);
fluctuation(:,1:res,3) = F_var;

for chan_num=1:N_chan
    [alpha(chan_num,1),rsq(chan_num,1), b_index(chan_num,1),conf(chan_num,1)] = routine_dfa_fit(n_vec, log10F_n(chan_num,:),res,1./F_var(chan_num,:));
end
fprintf('routine_dfa_calc> DFA exponents calculated. \n\n');
measures = [alpha,rsq,b_index,conf];

end






