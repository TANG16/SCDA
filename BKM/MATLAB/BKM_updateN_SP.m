function [N,  accept] = BKM_updateN_SP(fn_BKM_cov, N, theta, y, priorN)    
% Algorithm 2. Single-update Metropolis–Hastings algorithm 
% using the system process (i.e. binomial/Poisson) as the proposal distribution   

    T = size(N,2);
%     N1 = N(1,:);
%     Na = N(2,:);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    
    sigy = theta(9);
    
    accept = 0;
    % Prior update for the initial
    for t = 1:2    
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        N(1,t) = poissrnd(priorN(1)); % deltaN1 = 10
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t==1)
%             num = log(poisspdf(N(1,t),priorN(1)))...
%                  + log(binopdf(N(2,t+1),priorN(2),priorN(3)));
%             den = log(poisspdf(N1_old,priorN(1)))...
%                  + log(binopdf(N(2,t+1),priorN(2),priorN(3)));        
%             % Proposal is not symmetric:
%             num = num - log(poisspdf(N(1,t),priorN(1)));
%             den = den - log(poisspdf(N1_old,priorN(1)));
            % Or simply with cancelations:
%             num = log(binopdf(N(2,t+1),priorN(2),priorN(3)));
%             den = log(binopdf(N(2,t+1),priorN(2),priorN(3))); 
%% Double cancelation
            num = 1;
            den = 1;
        else
    %         num = log(poisspdf(N(1,t),priorN(1)))...
    %             + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
    %         den = log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)))...
    %             + log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t)));        
    %         % Proposal is not symmetric:
    %         num = num - log(poisspdf(N(1,t),priorN(1)));
    %         den = den - log(poisspdf(N1_old,priorN(1)));
            %% Or simply with cancelations:
            num = log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
            den = log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t))); 
        end
        
        % All other prior terms cancel in the acceptance probability.
        
        % Acceptance probability of MH step:
         A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
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
        
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t==1)   
%%    %         num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
    %              + log(binopdf(N(2,t),priorN(2),priorN(3))) ...
%                  + log(poisspdf(N(1,t+1),priorN(1)));
    %         den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
    %              + log(binopdf(Na_old,priorN(2),priorN(3))) ...
    %              + log(poisspdf(N(1,t+1),priorN(1)));         
    %         % Proposal is not symmetric:
    %         num = num - log(binopdf(N(2,t),priorN(2),priorN(3)));
    %         den = den - log(binopdf(Na_old,priorN(2),priorN(3)));
            % Or simply with cancelations:
%             num = - 0.5*((y(t)-N(2,t))^2)/sigy;
%             den = - 0.5*((y(t)-Na_old)^2)/sigy;
%% Or rather without y for t=1,2
              num = 1;
              den = 1;              
        else            
%%    %         num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
    %              + log(binopdf(N(2,t),priorN(2),priorN(3))) ...
    %              + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
    %         den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
    %              + log(binopdf(Na_old,priorN(2),priorN(3))) ...
    %              + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t)));         
    %         % Proposal is not symmetric:
    %         num = num - log(binopdf(N(2,t),priorN(2),priorN(3)));
    %         den = den - log(binopdf(Na_old,priorN(2),priorN(3)));
            % Or simply with cancelations:
%             num = - 0.5*((y(t)-N(2,t))^2)/sigy ...
%                  + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
%             den = - 0.5*((y(t)-Na_old)^2)/sigy ...
%                   + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t))); 
%% Or rather without y for t=1,2
            num = log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
            den = log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t))); 
        end
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));

        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
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
        
%%        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
%         num = log(poisspdf(N(1,t),N(2,t-1).*rho(t-1).*phi1(t-1)))...
%             + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
%         den = log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)))...
%             + log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t)));        
%         % Proposal is not symmetric:
%         num = num - log(poisspdf(N(1,t),N(2,t-1).*rho(t-1).*phi1(t-1)));
%         den = den - log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)));
        %% Or simply with cancelations:
        if (t==T)
            num = 1;
            den = 1;
        else
            num = log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
            den = log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t))); 
        end
        % All other prior terms  cancel in the acceptance probability.
   
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
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
        
%%        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
%         num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
%              + log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1))) ...
%              + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
%         den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
%              + log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1))) ...
%              + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t)));         
%         % Proposal is not symmetric:
%         num = num - log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1)));
%         den = den - log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1)));
        %% Or simply with cancelations:
        if (t==T)
            num = - 0.5*((y(t)-N(2,t))^2)/sigy;
            den = - 0.5*((y(t)-Na_old)^2)/sigy;
         else            
            num = - 0.5*((y(t)-N(2,t))^2)/sigy ...
                 + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
            den = - 0.5*((y(t)-Na_old)^2)/sigy ...
                  + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t)));   
        end
        % All other prior terms cancel in the acceptance probability.
       
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            N(2,t) = Na_old;
        end        
    end
end    