function [h, A_sum] = update_h(y, tau, h, theta, delta_h, hyper)

    T = length(y);
         
    h0 = theta(:,1);
%     g0 = theta(:,2);
    omegah = theta(:,3);
%     omegag = theta(:,4);

    accept = 0;
    A_sum = 0;
    
    for t = 1:T    
        % Keep a record of the current h value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t == 1)
            m = 0;
            s2 = hyper.Vh;
        else
            m = h(t-1);
            s2 = 1;
        end
        mean_y = tau(t);
        var_y = exp(h0 + omegah*h(t));
        var_y_old = exp(h0 + omegah*h_old);
   
%         num = -0.5*log(var_y) - 0.5*((y(t)-mean_y)^2)/var_y + ...
%               -0.5*log(s2) - 0.5*((h(t)-m)^2)/s2;
%         den = -0.5*log(var_y_old) - 0.5*((y(t)-mean_y)^2)/var_y_old + ...
%               -0.5*log(s2) - 0.5*((h_old-m)^2)/s2;
        num = -0.5*log(var_y) - 0.5*((y(t)-mean_y)^2)/var_y + ...
              - 0.5*((h(t)-m)^2)/s2;
        den = -0.5*log(var_y_old) - 0.5*((y(t)-mean_y)^2)/var_y_old + ...
              - 0.5*((h_old-m)^2)/s2;          

        if (t < T)   
            m_2 = h(t);
            m_2_old = h_old;

%             num = num - 0.5*(log(2*pi) + log(1) + ((h(t+1)-m_2)^2)/1);         
%             den = den - 0.5*(log(2*pi) + log(1) + ((h(t+1)-m_2_old)^2)/1);            
            num = num - 0.5*(((h(t+1)-m_2)^2)/1);         
            den = den - 0.5*(((h(t+1)-m_2_old)^2)/1);          
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
            h(t) = h_old;
        end        
    end
end
