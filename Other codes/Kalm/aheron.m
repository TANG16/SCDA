%mrr estimates for survival, from model phi1(Dec)/phia(Feb)/lambda(year)
Mu = [-0.5210;0.1559;0.2850;0.0189];

%variances-covariances
Sigma = [0.0097  0.0001  0.0051 -0.0018
         0.0001  0.0026  0.0002 -0.0003
         0.0051  0.0002  0.0217 -0.0019
        -0.0018 -0.0003 -0.0019  0.0042];

%number censused, yrs 1975-1990
Census = [5100;5075;5125;5375;5200;5575;5675;5350;5425;5700;5775;5109;5166;5637;5898;6143];


% covariates (optional)
%
%         Dec    Jan    Feb    Mar     Year of ringing = year of Dec Temp)
%         ^^^    ^^^    ^^^    ^^^     ^^^^
 Temp = [ 8.1    6.8    4.4    4.8    %1974
          5.3    5.9    4.5    4.8
          2.0    2.8    5.2    6.9
          6.1    3.4    2.8    6.7
          3.9   -0.4    1.2    4.7
          5.8    2.3    5.7    4.7
          5.6    4.9    3.0    7.9    %1980
          0.3    2.6    4.8    6.1
          4.4    6.7    1.7    6.4
          5.6    3.8    3.3    4.7
          5.2    0.8    2.1    4.7
          6.3    3.5   -1.1    4.9    %1985
          6.2    0.8    3.6    4.1
          5.6    5.3    4.9    6.4
          7.5    6.1    5.9    7.5
          4.9    6.5    7.3    8.3 ]; %1989


%census covariates
Decc = Temp(:,1); Janc = Temp(:,2); Febc = Temp(:,3); Marc = Temp(:,4);
Yrc = linspace(0,1,16)';

%scale Dec and Feb, as in ring-recovery exercise
Decc = Decc - 4.5909; %Dec - mean(Dec_1975:Dec_1985)
Febc = Febc - 2.8700; %Feb - mean(Feb_1976:Feb_1985)
