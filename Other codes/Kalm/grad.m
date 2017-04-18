function g = grad(funfcn,x,varargin)

%GRAD  Find gradient of function of a vector variable.
%   GRAD('F',X,P1,P2,...) approximates the gradient of F
%   at point X. The call provides for additional arguments
%   which are passed to the objective function, F(X,P1,P2,...).


% Convert to inline function as needed.
funfcn = fcnchk(funfcn,length(varargin));

%initialisations
delta = 10^(-6);
t = length(x);
g = zeros(t,1);
Dx = delta*eye(t);

%gradient
for i = 1:t
  g(i) = (feval(funfcn,x+Dx(:,i),varargin{:})-...
	 feval(funfcn,x-Dx(:,i),varargin{:}))/(2*delta);
end
