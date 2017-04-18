function h = hessian(funfcn,x,varargin)

%HESSIAN   Find hessian of function of a vector variable.
%   HESSIAN('F',X,P1,P2,...) approximates the hessian of F 
%   at point X. The call provides for additional arguments 
%   which are passed to the objective function, F(X,P1,P2,...). 


% Convert to inline function as needed.
funfcn = fcnchk(funfcn,length(varargin));

%initialisations
delta = 10^(-6);
t = length(x);
h = zeros(t);
Dx = delta*eye(t);

for i = 1:t
  for j = 1:t
    h(i,j)=(feval(funfcn,x+Dx(:,i)+Dx(:,j),varargin{:})-...
            feval(funfcn,x+Dx(:,i)-Dx(:,j),varargin{:})-...
            feval(funfcn,x-Dx(:,i)+Dx(:,j),varargin{:})+...
            feval(funfcn,x-Dx(:,i)-Dx(:,j),varargin{:}))/(4*delta^2);
  end
end
