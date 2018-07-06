function [tau, A_sum] = update_tau(y, tau, h, g, theta, delta_tau, hyper)

    T = length(y);
  
    h0 = theta(:,1);
    g0 = theta(:,2);
    omegah = theta(:,3);
    omegag = theta(:,4);

    accept = 0;
    A_sum = 0;
    
    for t = 1:T    
        % Keep a record of the current h value being updated
        tau_old = tau(t);
        % normal  RW for h
        tau(t) = tau_old + delta_tau*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t == 1)
            m = 0;
            s2 = hyper.Vtau*exp(g0 + omegag*g(t));
        else
            m = tau(t-1);
            s2 = exp(g0 + omegag*g(t));
        end
        mean_y = tau(t);
        mean_y_old = tau_old;
        var_y = exp(h0 + omegah*h(t));
    
%         num = -0.5*log(var_y) - 0.5*((y(t)-mean_y)^2)/var_y + ...
%               -0.5*log(s2) - 0.5*((h(t)-m)^2)/s2;
%         den = -0.5*log(var_y) - 0.5*((y(t)-mean_y_old)^2)/var_y + ...
%               -0.5*log(s2) - 0.5*((h_old-m)^2)/s2;
        num = - 0.5*((y(t)-mean_y)^2)/var_y + ...
              - 0.5*((tau(t)-m)^2)/s2;
        den = - 0.5*((y(t)-mean_y_old)^2)/var_y + ...
              - 0.5*((tau_old-m)^2)/s2;          

        if (t < T)   
            m_2 = tau(t);
            m_2_old = tau_old;
            var_tau = exp(g0 + omegag*g(t+1));
           
%             num = num - 0.5*(log(2*pi) + log(var_tau) + ((tau(t+1)-m_2)^2)/var_tau);         
%             den = den - 0.5*(log(2*pi) + log(var_tau) + ((tau(t+1)-m_2_old)^2)/var_tau);           
            num = num - 0.5*(((tau(t+1)-m_2)^2)/var_tau);         
            den = den - 0.5*(((tau(t+1)-m_2_old)^2)/var_tau);          
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
            tau(t) = tau_old;
        end        
    end
end
