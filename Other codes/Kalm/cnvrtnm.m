function a = cnvrtnm(x)
%global modelname X1name Xaname Xlamname X1 Xa Xlam phiatype
%global X1_C Xa_C Xlam_C Xa_A Xlam_T

% converts first character of string x to upper case
% and changes all commas, parentheses and colons
% to underscores

a = x;
if ~isstr(a), error('Argument must be a string'); end
if x(1)=='a', a(1) = 'A'; end
if x(1)=='c', a(1) = 'C'; end
if x(1)=='t', a(1) = 'T'; end
if x(1)=='v', a(1) = 'V'; end

n = length(x);

for i=1:n,
  if x(i)==',' | x(i)=='(' | x(i)==')' | x(i)==':',
    a(i) = '_';
  end
end
 
if a(n)=='_',
  a = a(1:n-1);
end
