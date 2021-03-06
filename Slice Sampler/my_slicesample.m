function [rnd, neval] = my_slicesample(initial,nsamples,varargin)
%SLICESAMPLE Slice sampling method.
%   RND = SLICESAMPLE(INITIAL,NSAMPLES,'pdf',PDF) draws NSAMPLES random
%   samples from a target distribution with the density function PDF using
%   the slice sampling. INITIAL is a row vector or scalar containing the
%   initial value of the random sample sequences. INITIAL must be within
%   the domain of the target distribution. NSAMPLES is the number of
%   samples to be generated. PDF is a function handle created using @. PDF
%   takes only one argument as an input and this argument has the same type
%   and size as INITIAL. It defines a function that is proportional to the
%   target density function. If log density function is preferred, 'pdf'
%   can be replaced with 'logpdf'. The log density function is not
%   necessarily normalized, either. 
%  
%   RND = SLICESAMPLE(...,'width',W) performs slice sampling for the target
%   distribution with a typical width W. W is a scalar or vector. If it is
%   a scalar, all dimensions are assumed to have the same typical widths.
%   If it is a vector, each element of the vector is the typical width of
%   the marginal target distribution in that dimension. The default value
%   of W is 10.
% 
%   RND = SLICESAMPLE(...,'burnin',K) generates random samples with values
%   between the starting point and the K-th point omitted in the generated
%   sequence, but keep points after that. K is a non-negative integer. The
%   default value of K is 0. 
%
%   RND = SLICESAMPLE(...,'thin',M) generates random samples with M-1 out of
%   M values omitted in the generated sequence. M is a positive integer.
%   The default value of M is 1.
%
%   [RND, NEVAL] = SLICESAMPLE(...) also returns NEVAL as the averaged
%   number of function evaluations occurred in the slice sampling. NEVAL is
%   a scalar.
%   
%   Example: 
%      % Define a function proportional to a multi-modal density
%      f = @(x) exp( -x.^2/2).*(1+(sin(3*x)).^2).* (1+(cos(5*x).^2));        
%      area = quad(f,-5,5);
% 
%      % Generate a sample based on this density
%      N = 2000;
%      x = slicesample(1,N,'pdf',f,'thin',5,'burnin',1000);  
% 
%      % Plot a histogram of the sample
%      [binheight,bincenter] = hist(x,50);
%      h = bar(bincenter,binheight,'hist');
%      set(h,'facecolor',[0.8 .8 1]);
% 
%      % Superimpose the f function scaled to have the same area
%      hold on 
%      xd = get(gca,'XLim');
%      xgrid = linspace(xd(1),xd(2),1000);
%      binwidth = (bincenter(2)-bincenter(1));
%      y = (N*binwidth/area) * f(xgrid);
%      plot(xgrid,y,'r','LineWidth',2)
%      hold off
%   
%   See also: MHSAMPLE, RAND, HIST, PLOT.

%   Slice sampling algorithm is summarized as 
%   (a) Starting form x0. 
%   (b) Draw a real value, y, uniformly from (0, f(x0)). 
%   (c) Find an interval, around x0, that contains all, or much, of the 
%       slice S ={x : f(x)>f(x0)}. 
%   (d) Draw the new point, x1, uniformly from the part of the slice within
%       this interval as a sample from this distribution. 
%   (e) Assign x1-->x0, and go back to step (b). 
%   (f) Repeat until the requested number of samples are obtained. 
%     
%   There are different ways of implementing step (c). We employ a strategy
%   called stepping-out and shrinking-in suggested by Neal(2003).
%
%   Reference: 
%     R. Neal (2003), Slice Sampling, Annals of Statistics, 31(3), p705-767
 
% Copyright 2005-2011 The MathWorks, Inc.

    narginchk(4,Inf);  % initial, nsamples and pdf or logpdf are required.


    initial = initial(:)';
    % parse the information in the name/value pairs 
    pnames = {'pdf' ,'logpdf', 'width','burnin','thin'};
    dflts =  {[] [] 10 0 1};
    [pdf,logpdf,width,burnin,thin] = ...
           internal.stats.parseArgs(pnames, dflts, varargin{:});

    % log density is used for numerical stability
    if isempty(logpdf)
        logpdf = @(x) mylog(pdf(x)); % log density
    end
    % MAXITER is used to limit the total number of iterations for each step.
    maxiter = 200;
    dim = size(initial,2); % dimension of the distribution
    rnd = zeros(nsamples,dim); % place holder for the random sample sequence
    neval = nsamples;  % one function evaluation is needed for each slice of the density.

    M = nsamples*thin + burnin;
    e = exprnd(1,M,1); % needed for the vertical position of the slice.
    % equivalent to u=rand(M,1); e = -log(u);

    RW = rand(M,dim); % factors of randomizing the width
    RD = rand(M,dim); % uniformly draw the point within the slice
    x0 = initial; % current value 

    % bool function indicating whether the point is inside the slice.
    inside = @(x,th) (logpdf(x) > th); 

    % update using stepping-out and shrinkage procedures.
    for i = 1-burnin:nsamples*thin   
        % A vertical level is drawn uniformly from (0,f(x0)) and used to define
        % the horizontal "slice".
        z = logpdf(x0) - e(i+burnin); % z = log(f(x0)*u), a random point under the graph
        % Will be the threshold in the function inside
        
        % An interval [xl, xr] of width w is randomly position around x0 and then
        % expanded in steps of size w until both size are outside the slice.   
        
        % The choice of w is usually tricky.  <<<<<<<<<<<<<<<<<<<<<<
        r = width.*RW(i+burnin,:); % random width/stepsize
        xl = x0 - r; 
        xr = xl + width; 
        iter = 0;

        % step out procedure is performed only when univariate samples are drawn.
        if dim==1 
            % step out to the left.
            while (inside(xl,z) && iter<maxiter)
                xl = xl - width;
                iter = iter +1;
            end       
            if (iter>=maxiter || any(xl<-sqrt(realmax))) % It takes too many iterations to step out.
                error(message('stats:slicesample:ErrStepout')) 
            end;
            neval = neval +iter;
            % step out to the right
            iter = 0;  
            while (inside(xr,z)) && iter<maxiter
                xr = xr + width;
                iter = iter+1;        
            end
            if iter>=maxiter || any(xr>sqrt(realmax)) % It takes too many iterations to step out.
                 error(message('stats:slicesample:ErrStepout')) 
            end
        end
        neval = neval + iter;

        % A new point is found by picking uniformly from the interval [xl, xr].
        xp = RD(i+burnin,:).*(xr-xl) + xl;

        % shrink the interval (or hyper-rectangle) if a point outside the
        % density is drawn.
        iter = 0;  
        while(~inside(xp,z))&& iter<maxiter  
            rshrink = (xp>x0);
            xr(rshrink) = xp(rshrink);
            lshrink = ~rshrink;
            xl(lshrink) = xp(lshrink);
            xp = rand(1,dim).*(xr-xl) + xl; % draw again
            iter = iter+1;
        end
        if iter>=maxiter % It takes too many iterations to shrink in.
            error(message('stats:slicesample:ErrShrinkin')) 
        end
        x0 = xp; % update the current value 
        
        if (i>0 && mod(i,thin)==0) % Assign value after the BurnIn and for every thin-th iteration
            rnd(i/thin,:) = x0;
        end
        neval = neval +iter;
    end

    neval = neval/(M); % averaged number of evaluations
end
%-------------------------------------------------
function  y = mylog(x)
% MYLOG function is defined to avoid Log of Zero warnings. 
y = -Inf(size(x));
y(x>0) = log(x(x>0));

 


