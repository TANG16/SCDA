clear all;

S = 10;
M = 10000;
BurnIn = 1000;
Y = zeros(S,36);

N1 = zeros(S,36);
Na = zeros(S,36);

NN_DA = zeros(36,M,S);
TH_DA = zeros(10000,9,S);

NN_HMM = zeros(36,M,S);
TH_HMM = zeros(10000,9,S);

for ii = 3:S
    s = RandStream('mt19937ar','Seed',ii);
    RandStream.setGlobalStream(s); 
    
    [N, y] = BKM_genarate;
    Y(ii,:) = y;
    N1(ii,:) = N(1,:);
    Na(ii,:) = N(2,:);  

    [NN_d, Theta_d] = BKM_try_DA_generate(M, BurnIn, false, ii);
    NN_DA(:,:,ii) = NN_d;
    TH_DA(:,:,ii) = Theta_d; 
    
    [NN_h, Theta_h] = BKM_try_HMM_generate(M, BurnIn, false, ii);
    NN_HMM(:,:,ii) = NN_h;
    TH_HMM(:,:,ii) = Theta_h; 
end


subplot(3,1,1)
plot(N1')
subplot(3,1,2)
plot(Na')
subplot(3,1,3)
plot(Y')