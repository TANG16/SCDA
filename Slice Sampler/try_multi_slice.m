function d = try_multi_slice(theta)
    % the first 12 dimensions can be specified by certain bounded uniform priors, 
    % and the last 3 dimensions can be specified by Beta(2,20).
    % x0 is the set of initial values
    % for the first 12 dimensions:
    % Posterior{i} = slicesample(x0(i), nsamples, 'pdf', @(x) unifpdf(x,0,200)); 
    % for the last 3 "Beta" dimensions, where i = 13, 14, 15
    % Posterior{i} = slicesample(x0(i), nsamples, 'pdf', @(x) betapdf(x,2,20));
 
    D = size(theta); %s hould be 15
    ok = (all(theta(1:12)<=200 & theta(1:12)>=0) & all(theta(13:15<=1 & theta(13:15)>=0))); 

    if (ok ~= 1)
        d = -Inf;
        return 
    end
   
    % log of uniform on [0,200] pdf(x) = I{[0,200]}/200 ==> -log(200)
    d = 12*(-log(200));
    d = d + sum(log(betapdf(theta(13:15),2,20)));    
end