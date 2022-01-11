function [IP,QP]=Forward_FDEM(conductivity,thick)
global freq
%
%% Sensor characteristics (S structure)
%
S.x=        5;                                                             % x-coordinate receiver (m)
S.y=        0;                                                             % y-coordinate receiver (m)
S.z=        0;                                                             % z-coordinate receiver (m) - positive z-axis pointed down
S.height=   0;                                                             % Height of transmitter (m)
S.mom=      3000;                                                          % Transmitter moment (A.m^2)
S.ori=      'ZZ';                                                          % Coil orientation (Transmitter(X,Y,Z),Receiver(X,Y,Z))
%
%% Model characteristics (M structure)
%
M.con = conductivity;                                                      % Conductivity of layer(s) (S/m)
M.thick = thick;                                                           % Layer(s) thickness (m)
layer = length(conductivity);                                              % Layer(s) (-)
M.sus = 4*pi*1e-7*ones(1,layer);                                           % Susceptibility of layer(s) (-)
M.perm = 8.85e-12*ones(1,layer);                                           % Permittivity of layer(s) (F/m)
for i=1:length(freq)
    S.freq = freq(i);                                                      % Frequency (Hz)
%
%% Calculate forward response (ppm)
%
    [FWD_IP,FWD_QP] = FDEM1DFWD_RC(S,M);                                   % Reflection coefficient approach
    IP(i) = FWD_IP;
    QP(i) = FWD_QP;
end
end