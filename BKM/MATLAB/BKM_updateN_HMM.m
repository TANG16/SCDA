function [N,  accept] = BKM_updateN_HMM(fn_BKM_cov, N, theta, y, deltaN, priorN, N_max)    
% Algorithm 1. Single-update Metropolis–Hastings algorithm 
% using a uniform proposal distribution   
    T = size(N,2);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    
    sigy = theta(9);
    
    accept = 0;
    for t = 1:2            
        % Keep a record of the current N value being updated
        Na_old = N(t);
        % Uniform RW for Na
        N(t) = Na_old + round(deltaN(1) * (2*rand-1)) ; %deltaN=10.5!
        if (N(t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            num = log(binopdf(N(t),priorN(1),priorN(2)));
            den = log(binopdf(Na_old,priorN(1),priorN(2)));
            % Proposal terms cancel since proposal distribution is symmetric.
            % All other prior terms cancel in the acceptance probability.

            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            N(t) = Na_old;
        end        
    end

    for t = 3:T    
        % Keep a record of the current N1 value being updated
        Na_old = N(t);
        % Uniform RW for Na
        N(t) = Na_old + round(deltaN(1) * (2*rand-1)); %deltaN=10.5!
        if (N(t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            num = BKM_loglik_N_HMM(y(t), N(t), N(t-1), N(t-2), phia(t-1), phi1(t-2), rho(t-2), sigy, N_max);
            den = BKM_loglik_N_HMM(y(t), Na_old, N(t-1), N(t-2), phia(t-1), phi1(t-2), rho(t-2), sigy, N_max);

            % Proposal terms cancel since proposal distribution is symmetric.
            % All other prior terms cancel in the acceptance probability. 
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end

        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            N(t) = Na_old;
        end        
    end
end    