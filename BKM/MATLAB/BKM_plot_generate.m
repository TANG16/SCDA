S = 10;
Y = zeros(S,36);

N1= zeros(S,36);
Na = zeros(S,36);


for ii = 1:S
    [N, y] = BKM_genarate;
    Y(ii,:) = y;
    N1(ii,:) = N(1,:);
    Na(ii,:) = N(2,:);   
end


subplot(3,1,1)
plot(N1')
subplot(3,1,2)
plot(Na')
subplot(3,1,3)
plot(Y')