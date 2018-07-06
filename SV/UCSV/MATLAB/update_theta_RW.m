    function [theta, accept] = update_theta_RW(y, tau, h, g, theta, hyper, delta) %, oldloglik)

    accept = zeros(1,4);

    % h0 ~ N(b0, Vh0);          Vh0 = 10.0;
    % g0 ~ N(c0, Vg0);          Vg0 = 10.0;
    % omegah ~ N(0,Vomegah);    Vomegah = 0.2
    % omegag ~ N(0,Vomegag);    Vomegag = 0.2

% hyper.h0 = [0, 10];
% hyper.g0 = [0, 10];
% hyper.omegah = 0.2;
% hyper.omegag = 0.2;
    
    % Cycle through each theta in turn 
    % and propose to update using random walk MH with uniform proposal density:
   
    for ii = 1:4
        % Keep a record of the current theta value being updated
        oldtheta = theta(ii);
        oldloglik = loglik_trend(y, tau, h, g, theta, hyper);
        
        % Propose a new value using a RW with uniform proposal density
        theta(ii) = theta(ii) + delta(ii)*randn; 
%         if ((ii == 1) || ((ii == 2) && (abs(theta(ii))<1)) || ((ii == 3) && (theta(ii)>0)))
            newloglik = loglik_trend(y, tau, h, g, theta, hyper);

            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:

            % Add in prior terms to the acceptance probability
            switch ii
                case 1
                    newprior = -0.5*((theta(ii)-hyper.h0(1)).^2)/hyper.h0(2);
                    oldprior = -0.5*((oldtheta-hyper.h0(1)).^2)/hyper.h0(2);
                case 2
                    newprior = -0.5*((theta(ii)-hyper.g0(1)).^2)/hyper.g0(2);
                    oldprior = -0.5*((oldtheta-hyper.g0(1)).^2)/hyper.g0(2);             
                case 3
                    newprior = -0.5*(theta(ii).^2)/hyper.omegah;
                    oldprior = -0.5*(oldtheta.^2)/hyper.omegah;                    
                otherwise    
                    newprior = -0.5*(theta(ii).^2)/hyper.omegag;
                    oldprior = -0.5*(oldtheta.^2)/hyper.omegag;                 
            end
            
            num = newloglik + newprior;
            den = oldloglik + oldprior;

            % All other prior terms (for other thetas) cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.

            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
%         else
%             A = 0;
%         end
        % To do the accept/reject step of the algorithm        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            accept(ii) = A;
        else  % Reject proposed move:
            % theta stays at current value:
            theta(ii) = oldtheta;
        end
    end  
    
end