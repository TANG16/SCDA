function y = ilogit(x)
% inverse logistic transform
y = (1+exp(-x)).^(-1);

