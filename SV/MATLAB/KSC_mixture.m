% Kim, Shephard and Chib (1998) RES, mixture approximation to the log standard
% normal density

KSC = [0.00730 -10.12999 5.79596
0.10556 -3.97281 2.61369
0.00002 -8.56686 5.17950
0.04395 2.77786 0.16735
0.34001 0.61942 0.64009
0.24566 1.79518 0.34023
0.25750 -1.08819 1.26261];

w_KSC = KSC(:,1);
m_KSC = KSC(:,2);
s2_KSC = SC(:,3);

const_KSC = -1.2704;

% y*_t = h_t + z_t
% y*_t = log(y_t^2+c)
% P(s_t = i| y*_t, h_t) ~ w_i * normpdf(y*_t | h_t + m_i - 1.2704, s2_i)