function [N,  accept] = BKM_updateN_HMM_bin_NEW(fn_BKM_cov, N, theta, y,...
                        deltaN, priorN, bin, logfact)    
% Algorithm 1. Single-update Metropolis–Hastings algorithm 
% using a uniform proposal distribution   
    T = size(N,2);
    
    [phi1, phia, rho] = fn_BKM_cov(theta);    
    sigy = theta(9);
    
%     accept = 0
    accept = zeros(1,T+8);
    for t = 1:2         
        % Keep a record of the current N value being updated
        Na_old = N(t);
        % Uniform RW for Na
%         N(t) = Na_old + round(deltaN(1) * (2*rand-1)) ; %deltaN=10.5!
        N(t) = binornd(priorN(2),priorN(3)); 
        if (N(t) > 0)
            %% Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
%             num = log(binopdf(N(t),priorN(1),priorN(2)));
%             den = log(binopdf(Na_old,priorN(1),priorN(2)));
%             num = N(t)*log(priorN(2)) + (priorN(1) - N(t))*log(1-priorN(2)) + ...
%                   logfact(priorN(1) + 1) - ...
%                   logfact(priorN(1) - N(t) + 1) - ...
%                   logfact(N(t) + 1);
%             den = Na_old*log(priorN(2)) + (priorN(1) - Na_old)*log(1-priorN(2)) + ...
%                   logfact(priorN(1) + 1) - ...
%                   logfact(priorN(1) - Na_old + 1) - ...
%                   logfact(Na_old + 1);      
            %% Or simply with cancelations:            
            num = N(t)*log(priorN(3)) - N(t)*log(1-priorN(3)) - ...
                  logfact(priorN(2) - N(t) + 1) - ...
                  logfact(N(t) + 1);
            den = Na_old*log(priorN(3)) - Na_old*log(1-priorN(3)) - ...
                  logfact(priorN(2) - Na_old + 1) - ...
                  logfact(Na_old + 1);       
              
              
%             if (t == 1)
% %                 loglam = log(priorN(1));  
%                 num = num + BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), priorN(1),...
%                       phia(t+1), 1, 1, sigy, bin, logfact);
%                 den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
%                       phia(t+1), phi1(t), rho(t), sigy, bin, logfact);  
%             end
              
            if (t == 2)  
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+1), N(t), priorN(1),...
                      phia(t), 1, 1, sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+1), Na_old, N(t-1),...
                      phia(t), 1, 1, sigy, bin, logfact);


                num = num + BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), N(t),...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);           
            end              
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
%             accept = accept+1;
            accept(t) = A;
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
                
            if (t == T)
                num = BKM_loglik_N_HMM_bin_v2(N(t), N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ... 
                      -  0.5*( ((y(t)-N(t)).^2)/sigy);
                den = BKM_loglik_N_HMM_bin_v2(Na_old, N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ...
                      -  0.5*(((y(t)-Na_old).^2)/sigy);   
                  
%                 num = num + BKM_loglik_N_HMM_bin_v2(N(t+1), N(t), N(t-1),...
%                       phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);
%                 den = den +  BKM_loglik_N_HMM_bin_v2(N(t+1), Na_old, N(t-1),...
%                       phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);
%                 num = num + BKM_loglik_N_HMM_EndCond(N(t-1), phi1(t-1), rho(t-1), bin, logfact);
%                 den = den + BKM_loglik_N_HMM_EndCond(N(t-1), phi1(t-1), rho(t-1), bin, logfact);

%                 num = num + BKM_loglik_N_HMM_bin_v2(1, 0, N(t),...
%                       phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
%                 den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
%                       phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
% %                 num = num + BKM_loglik_N_HMM_EndCond(N(t), phi1(t), rho(t), bin, logfact);
% %                 den = den + BKM_loglik_N_HMM_EndCond(Na_old, phi1(t), rho(t), bin, logfact);

            elseif (t == (T-1))
                num = BKM_loglik_N_HMM_bin_v2(N(t), N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ... 
                      -  0.5*( ((y(t)-N(t)).^2)/sigy);
                den = BKM_loglik_N_HMM_bin_v2(Na_old, N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ...
                      -  0.5*(((y(t)-Na_old).^2)/sigy); 
              
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+1), N(t), N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+1), Na_old, N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);    
                  
%                 num = num + BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), N(t),...
%                       phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
%                 den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
%                       phia(t+1), phi1(t), rho(t), sigy, bin, logfact);      
                num = num + BKM_loglik_N_HMM_bin_EndCond(N(t), phi1(t), rho(t), bin, logfact);
                den = den + BKM_loglik_N_HMM_bin_EndCond(Na_old, phi1(t), rho(t), bin, logfact);
            elseif (t == 3)
                num = BKM_loglik_N_HMM_bin_v2(N(t), N(t-1), priorN(1),...
                      phia(t-1), 1, 1, sigy, bin, logfact) ... 
                      -  0.5*( ((y(t)-N(t)).^2)/sigy);
                den = BKM_loglik_N_HMM_bin_v2(Na_old, N(t-1), priorN(1),...
                      phia(t-1), 1, 1, sigy, bin, logfact) ...
                      -  0.5*(((y(t)-Na_old).^2)/sigy); 
              
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+1), N(t), N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+1), Na_old, N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);    
                  
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), N(t),...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);                  

            else
                num = BKM_loglik_N_HMM_bin_v2(N(t), N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ... 
                      -  0.5*( ((y(t)-N(t)).^2)/sigy);
                den = BKM_loglik_N_HMM_bin_v2(Na_old, N(t-1), N(t-2),...
                      phia(t-1), phi1(t-2), rho(t-2), sigy, bin, logfact) ...
                      -  0.5*(((y(t)-Na_old).^2)/sigy); 
              
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+1), N(t), N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+1), Na_old, N(t-1),...
                      phia(t), phi1(t-1), rho(t-1), sigy, bin, logfact);    
                  
                num = num + BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), N(t),...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);
                den = den +  BKM_loglik_N_HMM_bin_v2(N(t+2), N(t+1), Na_old,...
                      phia(t+1), phi1(t), rho(t), sigy, bin, logfact);                  
            end
            

%                 num = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-N(2,t))^2)/(sigy)) ...
%                      + log(binopdf(N(2,t),N(1,t-1)+N(2,t-1),phia(t-1))) ...
%                      + log(binopdf(N(2,t+1),N(1,t)+N(2,t),phia(t))) ...
%                      + log(poisspdf(N(1,t+1),N(2,t).*rho(t).*phi1(t)));
% 
%                 den = - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na_old)^2)/(sigy)) ...
%                      + log(binopdf(Na_old,N(1,t-1)+N(2,t-1),phia(t-1))) ...
%                      + log(binopdf(N(2,t+1),N(1,t)+Na_old,phia(t))) ...
%                      + log(poisspdf(N(1,t+1),Na_old.*rho(t).*phi1(t))); 

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
%             accept = accept+1;
            accept(t) = A;
        else  % Reject proposed move:
            % Na stays at current value:
            N(t) = Na_old;
        end        
    end
end    