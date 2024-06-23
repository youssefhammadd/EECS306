clc
clear
close all
%parameter definition




Fs=1000;
dt=1/Fs;
t=0:dt:(1-dt);
fm=2;
Wm=2*pi*fm;
Am=4;
%modulating signal
m_t= Am .* cos(Wm *t);


figure();
plot(t,m_t);
xlabel('Time (msec)');
ylabel('m(t)');
title('Modulating Signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ac=2;
fc=10;
Wc=2*pi*fc;
%Unmodulating signal
c_t=Ac .* cos(Wc*t);
%generating the modulated signal
m_t_hat= Am .* sin(Wm*t);%quadrature of the Modulating signal
c_t_hat=Ac .*sin(Wc *t);%quadrature of the Unmodulating signal
%generating the USB and plot it
s_t_USB=m_t .*c_t - m_t_hat.* c_t_hat;
%generating the LSB 
s_t_LSB=m_t .*c_t + m_t_hat.* c_t_hat;

%plot of USB
figure();
plot(t,s_t_USB);
xlabel('Time (msec)');
ylabel('s(t) USB');
title('USB modulated signal');
%plot of LSB
figure();
plot(t,s_t_LSB);
xlabel('Time (msec)');
ylabel('s(t) LSB');
title('LSB modulated signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw the spectrum of both USB and LSB

f = -Fs/2:1:Fs/2 -1;  % Frequency vector

spectrum_LSB = fftshift(abs(fft(s_t_LSB/Fs))); % Shift zero frequency component to center
spectrum_USB = fftshift(abs(fft(s_t_USB/Fs))); % Shift zero frequency component to center

figure();
plot(f, spectrum_LSB);
xlim([-15 15]);
xlabel('Frequency (KHz)');
ylabel('Magnitude');
title('Spectrum of Modulated Signal LSB');
grid on;



figure();
plot(f, spectrum_USB);
xlim([-15 15]);
xlabel('Frequency (KHz)');
ylabel('Magnitude');
title('Spectrum of Modulated Signal USB');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%demodulator (coherent demodulator)

Demod_signal = s_t_LSB.*cos(Wc*t);
[b,a] = butter(4,0.01,'low');
Filtered_Demod_signal= filter(b,a,Demod_signal);
Spectrum_Filtered_Demod_signal = fftshift(fft(Filtered_Demod_signal));

figure();
plot(f,abs(Spectrum_Filtered_Demod_signal)/Fs);
title('Freq. Spectrum of Demodulated Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on
figure();
plot(t,Filtered_Demod_signal)
xlabel('Time(s)')
ylabel('Amplitude')
title('Demodulated Signal')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introducing a phase offset between carrier and local oscillator
phase_offset = pi / 4;  % Example phase offset of pi/4 radians
demodulated_signal_phase_offset = s_t_LSB .* cos(Wc * t + phase_offset);

% Demodulated signal after low-pass filtering with phase offset
demodulated_signal_filtered_phase_offset = filter(b, a, demodulated_signal_phase_offset);

% Plotting the demodulated signal with phase offset
figure();
plot(t, demodulated_signal_filtered_phase_offset);
title('Demodulated Signal with Phase Offset');
xlabel('Time (s)');
ylabel('Amplitude (V)');
grid on;