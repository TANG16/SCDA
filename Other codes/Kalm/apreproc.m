global Mu Sigma nmu Census ncensus

nmu = length(Mu);
ncensus = length(Census);

d = size(Sigma);
if any(d-nmu),
  fprintf('Sigma should be a square matrix of order %2.0f. ',nmu) 
  fprintf(['Check file ',datafilename,'.m\n\n']);
  break
elseif norm(Sigma-Sigma') > 1e-6
  fprintf(['Sigma should be a symmetric matrix. Check file ',datafilename,'.m\n\n']);
  break
else
  [r p] = chol(Sigma);
  if p > 0 
    fprintf(['Sigma should be a positive definite matrix. Check file ',datafilename,'.m\n\n']);
    break
  end
end
