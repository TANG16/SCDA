function [N, accept, A_sum] = BKM_updateN_U(fn_BKM_cov, N, theta, y, deltaN, priorN, logfact)    
% Algorithm 1. Single-update Metropolis–Hastings algorithm 
% using a uniform proposal distribution   
    
% While proposing N1(t) and Na(t) keep in mind the system restrictions:
% N1(t) + Na(t) >= Na(t+1)
% for N1(t)
% num = log(poisspdf(N(1,t),N(2,t-1).*rho(t-1).*phi1(t-1)))...
%       + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
% den = log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)))...
%       + log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t)));   
% hence:  N(1,t) + N(2,t) >= N(2,t+1)

% for Na(t)
% num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
%       + log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1))) ...
%       + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
% den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
% + log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1))) ...
% + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t))); 
% hence:  N(2,t)<= N(1,t-1) + N(2,t-1)

    T = size(N,2);
%     N1 = N(1,:);
%     Na = N(2,:);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    
    sigy = theta(9);
    
%     accept = 0;
    accept = zeros(1,T+T+8);  
    A_sum = 0;
    % RW update for the initial
    for t = 1:2    
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        N(1,t) = N1_old + round(deltaN(1) * (2*rand-1)) ; % deltaN1 = 10.5
        % Automatically reject any moves where N1<= 0
        if ((N(1,t) > 0) && (N(1,t) + N(2,t) >= N(2,t+1)))
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            % Add in prior (Poisson) terms to the acceptance probability
            if (t==1) 
                loglam = log(priorN(1));
%                 num =  - exp(loglam)  + N(1,t)*loglam - logfact(N(1,t) + 1);
%                 den =  - exp(loglam) + N1_old*loglam - logfact(N1_old + 1);   
                num =  N(1,t)*loglam - logfact(N(1,t) + 1);
                den =  N1_old*loglam - logfact(N1_old + 1);                   
    
%                 num = log(poisspdf(N(1,t),priorN(1)));
%                 den = log(poisspdf(N1_old,priorN(1)));     
            else
%                 num = log(poisspdf(N(1,t),priorN(1)))...
%                     + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
%                 den = log(poisspdf(N1_old,priorN(1))...
%                     + log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t))); 

                loglam = log(priorN(1));
                num = N(1,t)*loglam - logfact(N(1,t) + 1);
                den = N1_old*loglam - logfact(N1_old + 1);
                num = num + N(1,t)*log(1-phia(t)) + ...
                      logfact(N(1,t)+N(2,t) + 1) - ...
                      logfact(N(1,t)+N(2,t) - N(2,t+1) + 1);
                den = den + N1_old*log(1-phia(t)) + ...
                      logfact(N1_old+N(2,t) + 1) - ...
                      logfact(N1_old+N(2,t) - N(2,t+1) + 1); 
            end
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
        % Propose a new value for N1
        N(2,t) = Na_old + round(deltaN(2) * (2*rand-1)) ; % deltaNa = 100.5
        % Automatically reject any moves where N1<= 0
        if (N(2,t) > 0) 
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            if (t==1)
%                 num = log(binopdf(N(2,t),priorN(2),priorN(3)));
%                 den = log(binopdf(Na_old,priorN(2),priorN(3)));                
                num = N(2,t)*log(priorN(3)) -N(2,t)*log(1-priorN(3))  - ...
                      logfact(priorN(2) - N(2,t) + 1) - ...
                      logfact(N(2,t) + 1);                
                den = Na_old*log(priorN(3)) -Na_old*log(1-priorN(3)) - ...
                      logfact(priorN(2) - Na_old + 1) - ...
                      logfact(Na_old + 1);                 
            else
                if (N(2,t) <= N(1,t-1) + N(2,t-1))
    %%                 num = log(binopdf(N(2,t),priorN(2),priorN(3))) ...
    %                      + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
    %                 den = log(binopdf(Na_old,priorN(2),priorN(3))) ...
    %                      + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t)));         
    %                 loglam = log(N(2,t)) + log(rho(t)) + log(phi1(t));
    %                 num =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
    %                 num = num + N(2,t)*log(priorN(3)) + (priorN(2)-N(2,t))*log(1-priorN(3)) + ...
    %                       logfact(priorN(2) + 1) - ...
    %                       logfact(priorN(2) - N(2,t) + 1) - ...
    %                       logfact(N(2,t) + 1);
    %                 
    %                 loglam = log(Na_old) + log(rho(t)) + log(phi1(t));
    %                 den =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
    %                 num = num + Na_old*log(priorN(3)) + (priorN(2)-Na_old)*log(1-priorN(3)) + ...
    %                       logfact(priorN(2) + 1) - ...
    %                       logfact(priorN(2) - Na_old + 1) - ...
    %                       logfact(Na_old + 1);               
                    %% with cancelation
                    loglam = log(N(2,t)) + log(rho(t)) + log(phi1(t));
                    num =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
                    num = num + N(2,t)*log(priorN(3)) - N(2,t)*log(1-priorN(3))  - ...
                          logfact(priorN(2) - N(2,t) + 1) - ...
                          logfact(N(2,t) + 1);                
                    loglam = log(Na_old) + log(rho(t)) + log(phi1(t));
                    den =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
                    den = den + Na_old*log(priorN(3)) - Na_old*log(1-priorN(3)) - ...
                          logfact(priorN(2) - Na_old + 1) - ...
                          logfact(Na_old + 1); 
                else
                    num = -Inf;
                    den = Inf;
                end
            end
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

    for t = 3:T    
        %% N1 %%
        % Keep a record of the current N1 value being updated
        N1_old = N(1,t);
        % Propose a new value for N1
        N(1,t) = N1_old + round(deltaN(1) * (2*rand-1)) ; % deltaN1 = 10.5
        % Automatically reject any moves where N1<= 0
        if (N(1,t) > 0)
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            if (t==T)
%%                 num = log(poisspdf(N(1,t),N(2,t-1).*rho(t-1).*phi1(t-1)));
%                 den = log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)));                    
%                 loglam = log(N(2,t-1)) + log(rho(t-1)) + log(phi1(t-1));
%                 num =  - exp(loglam) + N(1,t)*loglam - logfact(N(1,t) + 1);
%                 den =  - exp(loglam) + N1_old*loglam - logfact(N1_old + 1);
                %% with cancelations
                loglam = log(N(2,t-1)) + log(rho(t-1)) + log(phi1(t-1));
                num = N(1,t)*loglam - logfact(N(1,t) + 1);                             
                den = N1_old*loglam - logfact(N1_old + 1);                
            else
                if ((N(1,t) + N(2,t) >= N(2,t+1)) && ((N1_old + N(2,t) >= N(2,t+1))))
    %%                 num = log(poisspdf(N(1,t),N(2,t-1).*rho(t-1).*phi1(t-1)))...
    %                     + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t)));
    %                 den = log(poisspdf(N1_old,N(2,t-1).*rho(t-1).*phi1(t-1)))...
    %                     + log(binopdf(N(2,t+1),N1_old+N(2,t),phia(t)));                             
    %                 loglam = log(N(2,t-1)) + log(rho(t-1)) + log(phi1(t-1));
    %                 num =  - exp(loglam) + N(1,t)*loglam - logfact(N(1,t) + 1);
    %                 den =  - exp(loglam) + N1_old*loglam - logfact(N1_old + 1);
    %                 num = num + N(1,t)*log(1-phia(t)) + ...
    %                       logfact(N(1,t)+N(2,t) + 1) - ...
    %                       logfact(N(1,t)+N(2,t) - N(2,t+1) + 1);
    %                 den = den + N1_old*log(1-phia(t)) + ...
    %                       logfact(N1_old+N(2,t) + 1) - ...
    %                       logfact(N1_old+N(2,t) - N(2,t+1) + 1);   
                    %% with cancelations
                    loglam = log(N(2,t-1)) + log(rho(t-1)) + log(phi1(t-1));
                    num = N(1,t)*loglam - logfact(N(1,t) + 1);
                    num = num + N(1,t)*log(1-phia(t)) + ...
                          logfact(N(1,t)+N(2,t) + 1) - ...
                          logfact(N(1,t)+N(2,t) - N(2,t+1) + 1);                
                    den = N1_old*loglam - logfact(N1_old + 1);
                    den = den + N1_old*log(1-phia(t)) + ...
                          logfact(N1_old+N(2,t) + 1) - ...
                          logfact(N1_old+N(2,t) - N(2,t+1) + 1);  
                else
                   num = -Inf;
                   den = Inf;
                end
            end     
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
        % Keep a record of the current N1 value being updated
        Na_old = N(2,t);
        % Propose a new value for N1
        N(2,t) = Na_old + round(deltaN(2) * (2*rand-1)) ; % deltaNa = 100.5
        % Automatically reject any moves where Na<= 0
        if ((N(2,t) > 0) && ((N(1,t-1) + N(2,t-1) >= N(2,t))) && ((N(1,t-1)+N(2,t-1)) >= Na_old))
            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            if (t==T)
%                 num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
%                      + log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1)));
%                 den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
%                      + log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1)));
                 
%                 loglam = log(N(2,t)) + log(rho(t)) + log(phi1(t));                 
                num =  - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy));
                num = num + N(2,t)*log(phia(t-1)) -N(2,t)*log(1-phia(t-1)) - ...
                      logfact(N(1,t-1)+N(2,t-1) - N(2,t) + 1) - ...
                      logfact(N(2,t) + 1);             
%                 loglam = log(Na_old) + log(rho(t)) + log(phi1(t));
                den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy));
                den = den + Na_old*log(phia(t-1)) - Na_old*log(1-phia(t-1)) - ...                      
                      logfact(N(1,t-1)+N(2,t-1) - Na_old + 1) - ...
                      logfact(Na_old + 1);                  
             else            
%%                 num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
%                      + log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1))) ...
%                      + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
%                 den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
%                      + log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1))) ...
%                      + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t)));                          
%                 loglam = log(N(2,t)) + log(rho(t)) + log(phi1(t));
%                 num =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
%                 num = num + N(2,t)*log(phia(t-1)) + (N(1,t-1)+N(2,t-1)-N(2,t))*log(1-phia(t-1)) + ...
%                       logfact(N(1,t-1)+N(2,t-1) + 1) - ...
%                       logfact(N(1,t-1)+N(2,t-1) - N(2,t) + 1) - ...
%                       logfact(N(2,t) + 1);
%                 
%                 loglam = log(Na_old) + log(rho(t)) + log(phi1(t));
%                 den =  - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
%                 den = den + Na_old*log(phia(t-1)) + (N(1,t-1)+N(2,t-1)-Na_old)*log(1-phia(t-1)) + ...
%                       logfact(N(1,t-1)+N(2,t-1) + 1) - ...
%                       logfact(N(1,t-1)+N(2,t-1) - Na_old + 1) - ...
%                       logfact(Na_old + 1);               
                 %% with cancelation
                loglam = log(N(2,t)) + log(rho(t)) + log(phi1(t));                 
                num =  - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) + ...
                      - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
                num = num + N(2,t)*log(phia(t-1)) -N(2,t)*log(1-phia(t-1)) - ...
                      logfact(N(1,t-1)+N(2,t-1) - N(2,t) + 1) - ...
                      logfact(N(2,t) + 1);             
                loglam = log(Na_old) + log(rho(t)) + log(phi1(t));
                den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) + ...
                      - exp(loglam) + N(1,t+1)*loglam - logfact(N(1,t+1) + 1);
                den = den + Na_old*log(phia(t-1)) - Na_old*log(1-phia(t-1)) - ...                      
                      logfact(N(1,t-1)+N(2,t-1) - Na_old + 1) - ...
                      logfact(Na_old + 1);                                         
            end
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