function D = distance_bte(Y)

% -------------------------------------------------------------------------
% Born to export
% 
% distance:
% - takes on parameters in Y.
% - calls 'main_bte' which solves the model and simulates it.
% - computes distance between moments implied by the model and those
%   estimated from the data.
%
% Written by Bernardo Diaz de Astarloa @ PSU February 2015.
% -------------------------------------------------------------------------

format long;

% Assign parameter values from Y.
ahf       = Y(1);              %#ok Success parameter, beta distribution
bhf       = Y(2);              %#ok Failure parameter, beta distribution
delta     = Y(3);              %#ok Exogenous match separation rate

% Call main routine (solution + simulation). main_lbd gives MM.
main_bte;

% Data and weights.
xdistData = [-4.386842; 1.524812];   % (m,v) of BGD exports sales distribution (lognormal fit).
haz1Data  = 0.563250;                % BGD match separation rate: 1 yr old (from CJ).
haz2Data  = 0.452657;                % BGD match separation rate: 2+ yr old (from CJ).

data = cat(1,xdistData,haz1Data,haz2Data);
W    = eye(size(data,1));            % Weight matrix.

% Compute loss function.
try
    error = data-MM;
    % D     = error'*W*error;
    D = norm(error)/norm(data);

    nanflag = isnan(D); 
    if nanflag>0;
        D = D*10; 
    end
        
    % Print Diagnostics
    compare = cat(2,data,MM);
    
    fprintf('==============================================================\n');
    fprintf('Loop preliminaries\n');
    fprintf('\n')
    fprintf('\nData - Model comparison:\n');
    disp(   num2str(compare));   
    fprintf('\nParameters values:\n');
    disp(   num2str(Y));
    fprintf('\nLoss (weighted):\n');
    disp(   num2str(D));
    %fprintf('\nLoss (unweighted):\n');
    %disp(   num2str(D_unw));
    fprintf('==============================================================\n');

catch ER 

    % Something went wrong: high loss
    disp('Something went wrong when computing metric: high loss.');
    D = 1e12;    
end 