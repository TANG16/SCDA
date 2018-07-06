function [g, A_sum] = update_g(tau, g, theta, delta_g, hyper)

    T = length(tau);
    
%     h0 = theta(:,1);
    g0 = theta(:,2);
%     omegah = theta(:,3);
    omegag = theta(:,4);

    accept = 0;
    A_sum = 0;
    
    for t = 1:T    
        % Keep a record of the current h value being updated
        g_old = g(t);
        % normal  RW for h
        g(t) = g_old + delta_g*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t == 1)
            m = 0;
            s2 = hyper.Vg;
            mean_tau = 0; 
            var_tau =  hyper.Vtau*exp(g0 + omegag*g(t));
            var_tau_old =  hyper.Vtau*exp(g0 + omegag*g_old);            
        else
            m = g(t-1);
            s2 = 1;
            mean_tau = tau(t-1);
            var_tau = exp(g0 + omegag*g(t));
            var_tau_old = exp(g0 + omegag*g_old);            
        end

%         num = -0.5*log(var_tau) - 0.5*((tau(t)-mean_tau)^2)/var_tau + ...
%               -0.5*log(s2) - 0.5*((g(t)-m)^2)/s2;
%         den = -0.5*log(var_tau_old) - 0.5*((tau(t)-mean_tau)^2)/var_tau_old + ...
%               -0.5*log(s2) - 0.5*((g_old-m)^2)/s2;
        num = -0.5*log(var_tau) - 0.5*((tau(t)-mean_tau)^2)/var_tau + ...
              -0.5*((g(t)-m)^2)/s2;
        den = -0.5*log(var_tau_old) - 0.5*((tau(t)-mean_tau)^2)/var_tau_old + ...
              -0.5*((g_old-m)^2)/s2;          

        if (t < T)   
            m_2 = g(t);
            m_2_old = g_old;
            
%             num = num - 0.5*(log(2*pi) + log(1) + ((g(t+1)-m_2)^2)/1);         
%             den = den - 0.5*(log(2*pi) + log(1) + ((g(t+1)-m_2_old)^2)/1);          
            num = num - 0.5*(((g(t+1)-m_2)^2)/1);         
            den = den - 0.5*(((g(t+1)-m_2_old)^2)/1);            
        end
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
            g(t) = g_old;
        end        
    end
end
