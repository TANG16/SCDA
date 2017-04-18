function x2 = chisq(Mobs,Mexp)
% Pearson chisquare value

[nr,nc] = size(Mobs);
x2 = 0;

for i=1:nr,
  for j=i:nc,
    x2 = x2 + ((Mobs(i,j)-Mexp(i,j))^2)/Mexp(i,j);
  end
end

