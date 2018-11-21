%% Calculate DFA exponent
% measures = [alpha,rsq,b_index]
% Sami Auno, 14.11.2018


function [measures,fluctuation] = routine_dfa_calc(DATA, minSampleWindow, maxSampleWindow, numTimeWindows, avgFlag)

% if minSampleWindow < 10
%     minSampleWindow = 10;
% end
[N_chan,N_data] = size(DATA);

% calcuale Hilbert transform -> envelope
fprintf('Calculating signal envelope... \n\n');
DATA_uenv = zeros(N_chan,N_data);
for i=1:N_chan
    [DATA_uenv(i,:),~] = envelope(DATA(i,:));
end
fprintf('Signal envelope calculated. \n\n');

% number of time windows in DFA
res = numTimeWindows;                   

% initializing and defining variables

n_vec = 2*floor(logspace(log10(minSampleWindow),log10(maxSampleWindow),res)/2);   % only even numbers, spread logarithmically
F_n = zeros(N_chan,res); % fluctuation as a function of n_vec
n_windows = zeros(N_chan,res);
alpha = zeros(N_chan,1);
rsq = zeros(N_chan,1);
b_index = zeros(N_chan,1);

% compute cumulative sum.
% DATA_cs = zeros(N_chan,N_data);
DATA_cs = cumsum(DATA_uenv,2);

fprintf('Estimating fluctuation function... \n\n');

% compute fluctuation for each n_vec
for j=1:res

    [F_n(:,j),n_windows(:,j)] = routine_fluctuation(DATA_cs, n_vec(j), N_data, N_chan, avgFlag);

end
fprintf('Fluctuation functions calculated. \n\n');

fprintf('Estimating DFA exponents... \n\n');
% linear regression analysis, compute alpha
fluctuation = zeros(N_chan,res,2);

log10F_n = log10(F_n);
n_vec = log10(n_vec);

fluctuation(:,:,1) = log10F_n;
fluctuation(:,:,2) = repmat(n_vec,[N_chan,1]);

for chan_num=1:N_chan
    [alpha(chan_num,1),rsq(chan_num,1), b_index(chan_num,1)] = routine_dfa_fit(n_vec, log10F_n(chan_num,:),res,n_windows(chan_num,:));
end
fprintf('DFA exponents calculated. \n\n');
measures = [alpha,rsq,b_index];

end