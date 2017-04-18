global Ringed Mobs unrec nrows ncols Census ncensus const Toe

[nrows ncols] = size(Mobs);
unrec = Ringed - sum(Mobs')';			  % unrecovered birds
ncensus = length(Census);
const = sum(lfact(Ringed)) - sum(sum(lfact(Mobs))) - sum(lfact(unrec));
Toe = toeplitz([0,ones(1,nrows-1)],zeros(1,ncols));

lik_max = ( const + sum(sum(Mobs(Mobs>0).*log(Mobs(Mobs>0)))) ...
            + sum(unrec(unrec>0).*log(unrec(unrec>0))) ...
            - sum(Ringed(Ringed>0).*log(Ringed(Ringed>0))) );
df = sum(sum(Mobs > 0));

% fprintf('Loglikelihood for maximal model is\n')
% fprintf(' lik_max = %10.4f on %4.0f d.f.\n\n',lik_max,df)
