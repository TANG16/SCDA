function [h, accept, A_sum] = update_h_HMM_shift(y, h, theta, delta_h,  bins, bin_midpoint, shift)
% if non-shitfed: differently than in the paper: integrate out the even h(t)'s and impute
% the odd ones
% if shitfed: as the paper: integrate out the odd h(t)'s and impute the even ones
    T = length(y);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    stdev_y = exp((mu+bin_midpoint)/2);  

    accept = 0;
    A_sum = 0;
    
    if ~shift
        index = 1:2:T; % integrate out h0 (initial) and impute h1
    else
        index = 2:2:T; % integrate out h1 
    end
    for t = index % 1:2:T       t=2:2:T     t=t+2
        %% RW IMPUTATION
        % Keep a record of the current N1 value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        num = SV_loglik_semi_HMM(y(t), h(t), theta, bins, bin_midpoint, true);
        den = SV_loglik_semi_HMM(y(t), h_old, theta, bins, bin_midpoint, true);
        % Proposal terms cancel since proposal distribution is symmetric.
        % All other prior terms cancel in the acceptance probability. 
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
 
        A_sum = A_sum + A;
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            h(t) = h_old;
        end              
    end
    
    if (odd && ~shift) %one more imputation if non-shifted; if shifted it has been already integrated out
        h_old = h(T+1);
        % normal  RW for h
        h(T+1) = h_old + delta_h*randn;
        if (h(T+1) > 0) % automatically reject negative volatilities
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            num = SV_loglik_semi_HMM(y(T+1), h(T+1), theta, bins, bin_midpoint, false);
            den = SV_loglik_semi_HMM(y(T+1), h_old, theta, bins, bin_midpoint, false);
            % Proposal terms cancel since proposal distribution is symmetric.
            % All other prior terms cancel in the acceptance probability. 
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end
        A_sum = A_sum + A;
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            h(T+1) = h_old;
        end        
    end
end
