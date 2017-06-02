tic
% G2 = zeros(N_max+1, T-2);
for t = 3:T
 %     P = zeros(N_max+1, t);

    G2(1:N_max,1) = exp(-exp(loglam(t-2)) + IND(1:N_max)'*loglam(t-2)...
        - logfact(IND(1:N_max)' + 1)); 
    G2(N_max+1,1) = max(0,1 - sum(G2(1:(N_max),t-2)));

%     IND_ok = IND((IND + Na(t-1) - Na(t)) > 0);
%     P(IND_ok+1) = exp(Na(t)*log(phia(t-1)) + (IND_ok + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
%                 logfact(IND_ok + Na(t-1) + 1) - ...
%                 logfact(IND_ok + Na(t-1) - Na(t) + 1) - ...
%                 logfact(Na(t) + 1));   
%     loglik = loglik + log(sum(G.*P)); % piecewise multiplication enough here
end
toc

tic
G = zeros(N_max+1, T-2);
IND = (0:N_max)';
G(1:N_max,:) = IND(1:N_max)*loglam(1:(T-2)); 
G(1:N_max,:) = bsxfun(@minus, G(1:N_max,:), exp(loglam(1:(T-2)))); 
G(1:N_max,:) = bsxfun(@minus, G(1:N_max,:), logfact(IND(1:N_max) + 1)'); 
G(1:N_max,:) = exp(G(1:N_max,:));
G(N_max+1,:) = max(0,1 - sum(G(1:(N_max),:),1));
toc