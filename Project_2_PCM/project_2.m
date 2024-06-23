clc
clear
close all
%parameter definition
Fs = 10000;
dt = 1/Fs;
Tmin=0;               %starting time
Tmax= 0.06666667;     %ending time
t = Tmin:dt:Tmax; %drawing two periods of the Message
frequency = 30;       %Frequency of the message in Hertz
W = 2*pi*frequency;   %Frequency of the message in rad/sec
Amplitude = 2;        %Amplitude of the message 
Number_of_Bits=input('Enter Number of bits per sample :  ');
Quantization_Levels = 2^Number_of_Bits; %number of quantization error

percentage_of_nyquist=input('Sampling rate to the Nyquist rate : ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the analog sinusoidal signal with frequency 30 HZ and Amplitude 2
%volts and plot it
Analog_Message= Amplitude .* sin(W*t); %sinusoidal signal
% Plot the analog signal
figure();
plot(t,Analog_Message);
xlabel('Time (sec)');
ylabel('Message');
title('Analog Signal with the samples');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sampling process
Nyquist_Sampling_Rate=2*frequency; %nyquist rate of the message


Sampling_Frequency= percentage_of_nyquist*Nyquist_Sampling_Rate;% sampling rate


Ts=1/Sampling_Frequency; %sampling time
nsMax=floor(Tmax/Ts);    %get the range of time samples
n=0:nsMax;
t_sampling=n*Ts;
Sampled_Message= Amplitude .* sin(W*t_sampling); %the sampled signal 


figure();
plot(t_sampling,Sampled_Message);
ylim([-2 2]);
xlabel('Time (sec)');
ylabel('Message');
title('Sampled Signal');
hold on
stem(t_sampling,Sampled_Message,"filled");% function to draw delta at the sampling time
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quantization process

Vmax = 2;
Vmin = -2;
Step_size = (Vmax - Vmin)/Quantization_Levels;
partition = Vmin:Step_size:Vmax; %Levels between Vmax and Vmin with different step size
                                 %here we took the levels as a mid-rise
                                 %Distinct endpoints of different ranges
codebook = Vmin-(Step_size/2):Step_size:Vmax+(Step_size/2);%Quantization value for each partition

[index , Quantized_value]=quantiz(Sampled_Message,partition,codebook); %Quantization process
%the index determines on which partition interval to each sample
%Quantized_value determines the output of the quantizer, which contains the quantization values of the input signal

Length1= length(index);           %the number of index are the number of samples
Length2= length(Quantized_value); %the number of Quantized_value are also the number of samples


for i=1:Length1
    if(index(i)~=0)  % To make index as binary decimal so started from 0 to N-1 levels
       index(i)=index(i)-1;
    end 
    i=i+1;
 end   
  for i=1:Length2                                           
     if(Quantized_value(i)==Vmin-(Step_size/2)) %To make quantize value in between the levels
         Quantized_value(i)=Vmin+(Step_size/2); %because it is a mid rise
     end
 end    

 figure();
 stem(Quantized_value,"filled");                       % Display the Quantized values
 title('Quantized Signal');
 ylabel('Amplitude');
 xlabel('Time');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %encoding process


 code=de2bi(index,'left-msb');      % Convert the decimal to binary of the indeces
 k=1;
for i=1:Length1
    for j=1:Number_of_Bits
        coded(k)=code(i,j);       % convert code matrix to a coded row vector
                                  % each row is a code for a sample
        j=j+1;                   
        k=k+1;
    end
    i=i+1;
end


number_of_plotted_bits = Length1*Number_of_Bits+2;

figure();
stairs(coded);                    % Display the encoded signal
grid on;
grid minor;
axis([0 number_of_plotted_bits -0.1 1.1]);  
title('Encoded Signal');
ylabel('Amplitude');
xlabel('Time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Demodulation process
qunt=reshape(coded,Number_of_Bits,length(coded)/Number_of_Bits);
index=bi2de(qunt','left-msb');         % Getback the index in decimal form
Quantized_value=Step_size*index+Vmin+(Step_size/2); % getback Quantized values

figure(); 
plot(Quantized_value);                               % Plot Demodulated signal
title('Demodulated Signal');
ylabel('Amplitude');
xlabel('Time');
%the function takes the sampling frquency and number of bits and the efficiency of the channel 
y= Bit_Rate_function(Sampling_Frequency,Number_of_Bits); 
disp('The bit rate :')
disp(y);


