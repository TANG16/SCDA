function dummy = prtupmat(M,d)

%File missing from eagle.tar - PTB (10.01.00)

[m n] = size(M);

if d == 0,
   for i = 1:m,
	for j = 1:i-1,
	 fprintf('    ')
	end
	for j = i:n,
	 fprintf('%4.0f',M(i,j))
	end
      fprintf('\n')
   end
elseif d == 2,
   for i = 1:m,
	for j = 1:i-1,
	 fprintf('	')
	end
	for j = i:n,
	 fprintf('%6.2f',M(i,j))
	end
      fprintf('\n')
   end
else
	for i = 1:m,
	for j = 1:n,
	 fprintf('%10.3f',M(i,j))
	end
	fprintf('\n')
   end
end
