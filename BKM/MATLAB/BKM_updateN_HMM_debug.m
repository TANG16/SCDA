function [N,  accept, A_sum] = BKM_updateN_HMM_debug(fn_BKM_cov, N, theta, y, deltaN, priorN, N_max, logfact)    
% Algorithm 1. Single-update Metropolis–Hastings algorithm 
% using a uniform proposal distribution   
    T = size(N,2);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    
%     sigy = theta(9);
    
%     accept = 0;
    accept = zeros(1,T+8);
    A_sum = 0;    
    for t = 1:T            
        % Keep a record of the current N value being updated
        Na_old = N(t);
        NN_old = N;
        % Uniform RW for Na
        N(t) = Na_old + round(deltaN(1) * (2*rand-1)) ; %deltaN=10.5!
        if (N(t) > 0)
            % Calculate the log(acceptance probability): 
            num =  BKM_statespace_HMM_debug(N, theta, y, ...
                phi1, phia, rho,...
                priorN, N_max, logfact);
            den =  BKM_statespace_HMM_debug(NN_old, theta, y, ...
                phi1, phia, rho,...
                priorN, N_max, logfact);               
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
%             accept = accept+1;
            accept(t) = A;
        else  % Reject proposed move:
            % Na stays at current value:
            N(t) = Na_old;
        end        
    end
end    