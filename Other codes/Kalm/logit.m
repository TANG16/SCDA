function y = logit(x)
% logistic transform
y = log(x./(1-x));
