function f = lfact(n)
sz = size(n);
for i = 1:sz(1)
  for j = 1:sz(2)
    f(i,j) = sum(log(1:n(i,j)));
  end
end
%f = sum(log(1:n));
