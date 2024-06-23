clc
clear 
close all

% Define parameters
sampling_rate = 1000;  % Sampling rate in hz
dt = 1/sampling_rate;
num_of_periods=2;      %number of periods of m(t)
t = 0:dt:(1-dt)*num_of_periods ;    %time interval

% Define the shift value
shift_value = 0.5;

% Calculate the shifted time axis
t_m = t + shift_value;

% Calculate the shifted sawtooth wave
modulating_signal = -2*(t_m - floor(t_m)) +1;


%parameters of the carrier
fc=10;                 %carrier frequency=10 KHz
Ac=1;                   %carrier amplitude
Wc=2*pi*fc;             %carrier frequency=2*pi*fc rad/sec

unmodulated_signal= Ac.*cos(Wc*t);
Ka=input('Enter the modulation index: ');
Am=Ka*Ac; %to get the new amplitude of the 
AM_modulated_signal=(Ac+Am*modulating_signal) .* cos(Wc*t) ;
                        %the DSB-LC modulated signal



%Plot the modulating signal waveform
figure();
plot(t, modulating_signal);
xlabel('Time (msec)');
ylabel('m(t)');
title('Modulating Signal');


 

figure();
plot(t, unmodulated_signal);
xlabel('Time (msec)');
ylabel('c(t)');
title('unmodulated Signal');


figure();
plot(t, AM_modulated_signal);
xlabel('Time (msec)');
ylabel('s(t)');
title('Modulated Signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generating an FM signal for the same carrier then Plot the FM signal

Kf=input('Enter the frequency deviation: ');
sampling_rate=100000;
FM_modulated_signal=fmmod(modulating_signal,fc,sampling_rate,Kf);

figure();
plot(t, FM_modulated_signal);
xlabel('Time (msec)');
ylabel('sfm(t)');
title('FM Modulated Signal');


