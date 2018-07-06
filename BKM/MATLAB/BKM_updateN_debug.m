function [N, accept, A_sum] = BKM_updateN_debug(f,stdT, N, theta, y, deltaN, priorN, logfact)    
% Algorithm 1. Single-update Metropolis–Hastings algorithm 
% using a uniform proposal distribution   
    
    T = size(N,2);
%     N1 = N(1,:);
%     Na = N(2,:);
    
%     [phi1, phia, rho] = fn_BKM_cov(theta);    
%     sigy = theta(9);
    
%     accept = 0;
    accept = zeros(1,T+T+8);  
    A_sum = 0;
    % RW update for the initial
    for t = 1:T
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        NN_old = N;
        % Propose a new value for N1
        N(1,t) = N1_old + round(deltaN(1) * (2*rand-1));
        % Automatically reject any moves where Na<= 0
        if (N(1,t) > 0) 
            num = BKM_statespace(N, theta, y, f, stdT, priorN, logfact);
            den = BKM_statespace(NN_old, theta, y, f, stdT, priorN, logfact);
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
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
            % N1 stays at current value:
            N(1,t) = N1_old;
        end
        
        %% Na %%
        % Keep a record of the current Na value being updated
        Na_old = N(2,t);
        NN_old = N;
        % Propose a new value for Na
        N(2,t) = Na_old + round(deltaN(2) * (2*rand-1)) ; % deltaNa = 100.5
        % Automatically reject any moves where Na<= 0
        if (N(2,t) > 0) 
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            num = BKM_statespace(N, theta, y, f, stdT, priorN, logfact);
            den = BKM_statespace(NN_old, theta, y, f, stdT, priorN, logfact);
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
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
            accept(T+t) = A;                        
        else  % Reject proposed move:
            % N1 stays at current value:
            N(2,t) = Na_old;
        end        
    end
end    