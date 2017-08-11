function [h, accept, A_sum] = update_h_HMM(y, h, theta, delta_h,  bins, bin_midpoint)
%vintegrate out the odd h(t)'s and impute the even ones
    T = length(y);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    stdev_y = exp((mu+bin_midpoint)/2);  

    accept = 0;
    A_sum = 0;
    
    h0 = mu; %mu/(1-phi); % unconditional mean
    
    for t = 2:2:T % 1:2:T    
        %% RW IMPUTATION
        % Keep a record of the current N1 value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t == 2)
            Y = [y(t-1),y(t),y(t+1)];
            H = [h0,h(t),h(t+2)];
            H_old = [h0,h_old,h(t+2)];
         elseif ((t == T) && ~odd) % the last even y is indeed the last one
            Y = [y(t-1),y(t),NaN];
            H = [h(t-2),h(t),NaN];
            H_old = [h(t-2),h_old,NaN];   
         elseif ((t == T) && odd) % the last even y is not the last one, there is still the really last odd y
                                 % but there is no more h(t+2) to impute
            Y = [y(t-1),y(t),y(t+1)];
            H = [h(t-2),h(t),NaN];
            H_old = [h(t-2),h_old,NaN];   
         else
            Y = [y(t-1),y(t),y(t+1)];
            H = [h(t-2),h(t),h(t+2)];
            H_old = [h(t-2),h_old,h(t+2)]; 
         end
        num = SV_loglik_semi_HMM_v2(Y, H, theta, bins, bin_midpoint);
        den = SV_loglik_semi_HMM_v2(Y, H_old, theta, bins, bin_midpoint);
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
