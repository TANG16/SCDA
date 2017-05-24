function [N, newlikhood, accept] = BKM_updateN_v4(fn_BKM_calclikhood,  fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN)    
    
    T = size(N,2);
%     N1 = N(1,:);
%     Na = N(2,:);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    

    accept = 0;
    % RW update for the initials
    for t = 1:2    
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
%         N(1,t) = N1_old + round(deltaN(1) * (2*rand-1)) ; % deltaN1 = 10
        N(1,t) = poissrnd(priorN(1)); % deltaN1 = 10
        % Automatically reject any moves where N1<= 0
%         if (N(1,t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            newlikhood = fn_BKM_calclikhood(N);
            % Proposal distribution is not symmetric so  need to account for it        
            num = newlikhood - log(poisspdf(N(1,t),priorN(1)));
            den = oldlikhood - log(poisspdf(N1_old,priorN(1)));
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
%         else 
%             A = 0;
%         end
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
        
        %% Na %%
        % Keep a record of the current N1 value being updated
        Na_old = N(2,t);
        % Propose a new value for N1
%         N(2,t) = Na_old + round(deltaN(2) * (2*rand-1)) ; % deltaNa = 100
        N(2,t) = binornd(priorN(2),priorN(3)); % deltaNa = 100
        % Automatically reject any moves where N1<= 0
%         if (N(2,t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            newlikhood = fn_BKM_calclikhood(N);
            % Prior (Binomial) terms are already in the state space
            % likelhood
            % Proposal distribution is not symmetric so  need to account for it        
            num = newlikhood - log(binopdf(N(2,t),priorN(2),priorN(3)));
            den = oldlikhood - log(binopdf(Na_old,priorN(2),priorN(3)));
            % All other prior terms cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
%         else 
%             A = 0;
%         end
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
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        lam = N(2,t-1).*rho(t-1).*phi1(t-1);
        N(1,t) = poissrnd(lam);
        
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = fn_BKM_calclikhood(N);
        
        % Usually: add in prior (Poisson) terms to the acceptance
        % probability --> here already in the likelihood
        % Proposal distribution is not symmetric so  need to account for it
        num = newlikhood - log(poisspdf(N(1,t),lam));
        den = oldlikhood - log(poisspdf(N1_old,lam));
        % All other prior terms  cancel in the acceptance probability.
   
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
        
        %% Na %%
        % Keep a record of the current N1 value being updated
        Na_old = N(2,t);
        % Propose a new value for N1
        N(2,t) = binornd(N(1,t-1)+N(2,t-1),phia(t-1));
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = fn_BKM_calclikhood(N);
   
        % Usually: add in prior (Poisson) terms to the acceptance
        % probability --> here already in the likelihood
        num = newlikhood - log(binopdf(N(2,t),priorN(2),priorN(3)));
        den = oldlikhood - log(binopdf(Na_old,priorN(2),priorN(3)));
        % All other prior terms cancel in the acceptance probability.
       
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            N(2,t) = Na_old;
        end        
    end
end    