clear all;

arg1 = 1; % 1; % DA Kim
arg2 = 0; % 2; % DA RW eff
arg3 = 0; % semi DA
arg4 = 0; % semi DA shifting
arg5 = 0; % semi DA eff
arg6 = 1; % semi DA adaptive

% for ii = 1:6
ii=0; % 0: simulation; 1-6 different data sets
    try_SV_param_linux(ii,arg1,arg2,arg3,arg4,arg5,arg6)
% end    

% 
% arg0 = 4;
% arg1 = 1;
% arg2 = 2;
% arg3 = 0;
% arg4 = 0;
% arg5 = 1;
% arg0 = 1;

try_SV_param_linux(0,0,0,0,1,0)