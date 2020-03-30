%% Routine function for calculating DFA exponent alpha
% and other variables.
% [alpha,Rsqrd, b_index,conf] = routine_dfa_fit(n_vec,F_n,res,weights)


% Sami Auno, 2020

function [alpha,Rsqrd, b_index,conf] = routine_dfa_fit(n_vec,F_n,res,weights)
    
    % calculate linear and parabolic regression
%     reg1 = polyfit(n_vec,F_n,1);
    X = [ones(size(n_vec')) n_vec'];
    [reg1,stdX] = lscov(X,F_n',weights');
    reg1 = flipud(reg1)';
    reg2 = polyfit(n_vec,F_n,2);

    % calculare parameters for R-squared
    Ffit1 = polyval(reg1,n_vec);
    Ffit2 = polyval(reg2,n_vec);

    Fresid = F_n-Ffit1;
    SSresid = sum(Fresid.^2);
    SStotal = (double(res)-1) * var(F_n);

    % calculate parameters for parabolic index b
    E1 = rms(Ffit1-F_n);
    E2 = rms(Ffit2-F_n);
    
    conf = 1.96*stdX(2);
    Rsqrd = 1 - SSresid/SStotal;
    b_index = 1 - E2/E1;
    alpha = reg1(1);
end