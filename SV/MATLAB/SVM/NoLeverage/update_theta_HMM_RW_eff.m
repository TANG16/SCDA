function [theta, accept, oldloglik] = update_theta_HMM_RW_eff(y, h, theta,...
    bins, bin_midpoint, hyper, delta, oldloglik)

    accept = zeros(1,4); 

    % Cycle through each theta in turn 
    % and propose to update using random walk MH with uniform proposal density:
   
    for ii = 1:3
        % Keep a record of the current theta value being updated
        oldtheta = theta(ii);
%         oldloglik = loglik_h_HMM(y, h, theta, bins, bin_midpoint);
        
        % Propose a new value using a RW with uniform proposal density
        theta(ii) = theta(ii) + delta(ii)*randn; 
        if ((ii == 1) || ((ii == 2) && (abs(theta(ii))<1)) || ((ii == 3) && (theta(ii)>0)))
            newloglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint);

            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:

            % Add in prior terms to the acceptance probability
            switch ii
                case 1
                    newprior = -0.5*(theta(ii).^2)/hyper.M;
                    oldprior = -0.5*(oldtheta.^2)/hyper.M;                 
                case 2
                    newprior = (hyper.P(1)-1)*log((1+theta(ii))/2) + (hyper.P(2)-1)*log((1-theta(ii))/2); 
                    oldprior = (hyper.P(1)-1)*log((1+oldtheta)/2) + (hyper.P(2)-1)*log((1-oldtheta)/2);                 
                otherwise    
                    newprior = - (hyper.S(1)+1)*log(theta(ii)) - hyper.S(2)./theta(ii);
                    oldprior = - (hyper.S(1)+1)*log(oldtheta) - hyper.S(2)./oldtheta;                
            end
            num = sum(newloglik) + newprior;
            den = sum(oldloglik) + oldprior;

            % All other prior terms (for other thetas) cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.

            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else
            A = 0;
        end
        % To do the accept/reject step of the algorithm        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldloglik = newloglik;          
            accept(ii) = A;
        else  % Reject proposed move:
            % theta stays at current value:
            theta(ii) = oldtheta;
        end
    end  
    
    %%
    ii = 4;
     % Keep a record of the current theta value being updated
    oldtheta = theta(ii);
    oldloglik_y = loglik_y_HMM_eff(y, h, theta, bins, bin_midpoint);

    % Propose a new value using a RW with uniform proposal density
    theta(ii) = theta(ii) + delta(ii)*randn; 
    newloglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint);
    newloglik_y = loglik_y_HMM_eff(y, h, theta, bins, bin_midpoint);

    newprior = -0.5*(theta(ii).^2)/hyper.B;
    oldprior = -0.5*(oldtheta.^2)/hyper.B;      

    num = sum(newloglik) + newloglik_y + newprior;
    den = sum(oldloglik) + oldloglik_y + oldprior;

    % All other prior terms (for other thetas) cancel in the acceptance probability.
    % Proposal terms cancel since proposal distribution is symmetric.

    % Acceptance probability of MH step:
    A = min(1,exp(num-den));

    % To do the accept/reject step of the algorithm        
    % Accept the move with probability A:
    if (rand <= A)  % Accept the proposed move:
        % Update the log(likelihood) value:
        accept(ii) = A;
        oldloglik = newloglik;                  
    else  % Reject proposed move:
        % theta stays at current value:
        theta(ii) = oldtheta;
    end
end