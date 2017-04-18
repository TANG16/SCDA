function dummy = prtupsgn(M)

%File missing from eagle.tar - PTB (10.01.00)

[m n] = size(M);

for i = 1:m,
  for j = 1:i-1,
	fprintf('  ')
  end
	for j = i:n,
	if M(i,j)==1,
	fprintf('+ ')
	elseif M(i,j)==-1,
	fprintf('- ')
    else
	fprintf('. ')
    end
  end
	fprintf('\n')
end
