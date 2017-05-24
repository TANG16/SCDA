function [N, newlikhood, accept] = BKM_updateN_v2(fn_BKM_calclikhood,  fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN)    
    
    T = size(N,2);
%     N1 = N(1,:);
%     Na = N(2,:);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    

    accept = 0;
    
    %% N1 %%
    % RW update for the initial
    for t = 1:2    
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        N(1,t) = N1_old + round(deltaN(1) * (2*rand-1)) ; % deltaN1 = 10
        % Automatically reject any moves where N1<= 0
        if (N(1,t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            newlikhood = fn_BKM_calclikhood(N);
            % Add in prior (Poisson) terms to the acceptance probability
            num = newlikhood + log(poisspdf(N(1,t),priorN(1)));
            den = oldlikhood + log(poisspdf(N1_old,priorN(1)));
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % N1 stays at current value:
            N(1,t) = N1_old;
        end
    end

    for t = 3:T    
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        lam = N(2,t-1).*rho(t-1).*phi1(t-1);
        N(1,t) = poissrnd(lam);
        
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = fn_BKM_calclikhood(N);
        
        % Usually: add in prior (Poisson) terms to the acceptance probability
        % Proposal distribution is not symmetric so usually need to account for it
        % But here it is equal to the prior (the system process) 
        % so we even do not need to add prior/transition terms
        % [Keep in case the proposal changes]
        num = newlikhood; % + log(poisspdf(N1(t),lam));
        den = oldlikhood; % + log(poisspdf(N1_old,lam));
        % All other prior terms  cancel in the acceptance probability.
        % The acceptance ratio simplifies to the likelihood ratio! 
        num = num; % + log(poisspdf(N1_old,lam));
        den = den; % + log(poisspdf(N1(t),lam));
        
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % N1 stays at current value:
            N(1,t) = N1_old;
        end   
    end
    
    %% Na %%
    for t = 1:2    
        % Keep a record of the current N1 value being updated
        Na_old = N(2,t);
        % Propose a new value for N1
        N(2,t) = Na_old + round(deltaN(2) * (2*rand-1)) ; % deltaNa = 100
        % Automatically reject any moves where N1<= 0
        if (N(2,t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            newlikhood = fn_BKM_calclikhood(N);
            % Add in prior (Binomial) terms to the acceptance probability
            num = newlikhood + log(binopdf(N(2,t),priorN(2),priorN(3)));
            den = oldlikhood + log(binopdf(Na_old,priorN(2),priorN(3)));
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % N1 stays at current value:
            N(2,t) = Na_old;
        end        
    end

    for t = 3:T           
        % Keep a record of the current N1 value being updated
        Na_old = N(2,t);
        % Propose a new value for N1
        N(2,t) = binornd(N(1,t-1)+N(2,t-1),phia(t-1));
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = fn_BKM_calclikhood(N);
   
        % Usually: add in prior (Poisson) terms to the acceptance probability
        % Proposal distribution is not symmetric so usually need to account for it
        % But here it is equal to the prior (the system process) 
        % so we even do not need to add prior/transition terms
        % [Keep in case the proposal changes]
        num = newlikhood; % + log(binopdf(Na(t),priorN(2),priorN(3)));
        den = oldlikhood; % + log(binopdf(Na_old,priorN(2),priorN(3)));
        % All other prior terms  cancel in the acceptance probability.
        % The acceptance ratio simplifies to the likelihood ratio!        
        num = num; % + log(binopdf(Na_old,priorN(2),priorN(3)));
        den = den; % + log(binopdf(Na(t),priorN(2),priorN(3)));

        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % N1 stays at current value:
            N(2,t) = Na_old;
        end        
    end
end    