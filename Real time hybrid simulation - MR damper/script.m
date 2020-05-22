close all
clear all
clc
% This script runs the input parameters necessary for the MR Damper. A
% simulink file is included with this script, named "MRsim".
%% Simulink Parameters
    dt = 1/10000; % sec
    duration = 100;
    X = bodeoptions;
    X.FreqUnits = 'Hz';
    X.Grid = 'on';
    X.Xlim = [.1,10];
%% Model Parameters for MR Damper - Bouc-Wen Model
    A = 300;
    n = 2.0;
    c0a = 0.080; % kN.sec/mm
    c0b = 0.32; % kN.sec/mm
    c0c = 1.5/A; % kN.sec/mm
    k0 = 0.0; % kN/mm
    k1 = 0.0; % kN/mm
    x0 = 0.0; % mm
    c1a = 3.0; % kN.sec/mm
    c1b = 15.0; % kN.sec/mm
    c1c = 2.0/A; % kN.sec/mm
    alphaa = 0.11; % kN/mm
    alphab = 0.55; % kN/mm
    alphac = 1.0/A; % kN/mm
    betaa = 0.05; % mm^-2
    betab = 0.0020; % mm^-2
    betac = 5.2/A; % mm^-2

%% User specified Parameters
    % desired current: %%%%% Adjustable Parameter%%%%%
    id = 1*A; % Amps
    % effective current;
    ic = id;

%% Current dependent variables:
    Alpha = alphab + (alphaa - alphab)*exp(-alphac*ic);
    C0 = c0b + (c0a - c0b)*exp(-c0c*ic);
    C1 = c1b + (c1a - c1b)*exp(-c1c*ic);
    Beta = betab + (betaa - betab)*exp(-betac*ic);
    Gamma = Beta;

%% Numerical model for frame building
    M = 98300*eye(3); % kg
    K = [12.0 -6.84  0;
         -6.84 13.7 -6.84;
         0     -6.84 6.84]*10^7; % N/m
    [V,E] = eig(K,M);
    E = sqrt(E)/2/pi;
    Mr = V'*M*V;
    Kr = V'*K*V;
    z = .01;
    Cr = zeros(3,3);
    Cr(1,1) = 2*Mr(1,1)*z*E(1,1);
    Cr(2,2) = 2*Mr(2,2)*z*E(2,2);
    Cr(3,3) = 2*Mr(3,3)*z*E(3,3);
    C = inv(V')*Cr*inv(V); % Ns/m
    G = [-1 0 0]';
    I = [1 1 1]';

%% State-space model for frame building
    A1 = [zeros(3,3) eye(3);
         -M\K -M\C];
    B1 = [zeros(3,1);
          inv(M)*G];
    B2 = [zeros(3,1);
          -I];
    B = [B2 B1];
    C1 = [eye(3)*1000 zeros(3,3);
          -M\K -M\C];
    D = [B1 zeros(6,1)];
    bdg = ss(A1,B,C1,D);
    % bode(bdg,X)
clear B1 B2 A1 B C1 D

%% Load Ground Motion
    elcentro = load('elcentro.mat');
    kobe = load('kobe.mat');
    northridge = load('northridge.mat');
    elcentro = elcentro.elcentro;
    kobe = kobe.acc;
    northridge = northridge.acc;
    % choose ground motion: 
    gmotion = elcentro;
    duration = gmotion(end,1);
%% Actuator Model
    s = tf('s');
    Gact = 1.6e7/((s+151.7)*(s^2+250.4*s+1.061e5));
%% Lowpass filter
    fc = 5;
    [b,a] = butter(10,2*pi*fc,'s');
    lp = tf(b,a);
%% Run Simulation
    amp = 2;
    fsin = 10;
    sim('MRsim')
%% Figures,
% Displacement
    figure,
        subplot(3,1,1)
            plot(time, dout_RTHS(:,1), time, dout(:,1))
            grid on
            xlabel('Time (s)'); ylabel('Disp. (mm)')
            legend('w/ damper','no damper')
            title('First Floor')
            
        subplot(3,1,2)
            plot(time, dout_RTHS(:,2), time, dout(:,2))
            grid on
            xlabel('Time (s)'); ylabel('Disp. (mm)')
            title('Second Floor')
        subplot(3,1,3)
            plot(time, dout_RTHS(:,3), time, dout(:,3))
            grid on
            xlabel('Time (s)'); ylabel('Disp. (mm)')
            title('Third Floor')
% Acceleration
    figure,
        subplot(3,1,1)
            plot(time, aout_RTHS(:,1), time, aout(:,1))
            grid on
            xlabel('Time (s)'); ylabel('Acc. (m/s^{2})')
            legend('w/ damper','no damper')
            title('First Floor')
        subplot(3,1,2)
            plot(time, aout_RTHS(:,2), time, aout(:,2))
            grid on
            xlabel('Time (s)'); ylabel('Acc. (m/s^{2})')
            title('Second Floor')
        subplot(3,1,3)
            plot(time, aout_RTHS(:,3), time, aout(:,3))
            grid on
            xlabel('Time (s)'); ylabel('Acc. (m/s^{2})')
            title('Third Floor')
            
% Restoring Force hysteresis
    figure,
        plot(dout_RTHS(:,1), Fout_RTHS)
            grid on
            grid minor
            xlabel('1st Floor Disp. (mm)'); ylabel('Restoring Force (kN)')
            title('Restoring Force Hysteresis')

