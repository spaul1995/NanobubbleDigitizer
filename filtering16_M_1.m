clear all
Xk = zeros(1000000,10);
X2k = zeros(1000000,10);
for i = 1:10
filename = sprintf('data_%i.csv',fix(i));
% filename = 'data_1.csv';
start = 1;
ending = 1000000;
X1 = csvread(filename);

Xk(:,i)=X1(:,5);
X2k(:,i)=X1(:,11);
end

X = cat(1,Xk(:,1),Xk(:,2),Xk(:,3),Xk(:,4),Xk(:,5),Xk(:,6),Xk(:,7),Xk(:,8),Xk(:,9),Xk(:,10));
X2 = cat(1,X2k(:,1),X2k(:,2),X2k(:,3),X2k(:,4),X2k(:,5),X2k(:,6),X2k(:,7),X2k(:,8),X2k(:,9),X2k(:,10));




downsample_fac = 1;

Fs = 2500E6;                                       % Sampling frequency
T = 1/Fs;                                           % Sampling period
L = (10*ending-start)+1;                                           % Length of signal
t = (0:L-1)*T*1E6;                                      % Time vector
t=t';
Fn = Fs/2;                                          % Nyquist Frequency

X = downsample(X,downsample_fac);
X2 = downsample(X2,downsample_fac);
t = downsample(t,downsample_fac);
FX = fft(X)/L; 
FX2 = fft(X2)/L;% Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
Iv = 1:length(Fv);                                  % Index Vector
% figure(1)
% plot(Fv, abs(FX(Iv))*2)
% grid
% title('Fourier Transform Of Original Signal ‘X’')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
FXdcoc2 = fft(X2-mean(X2))/L; 
FXdcoc = fft(X-mean(X))/L;                          % Fourier Transform (D-C Offset Corrected)
% figure(2)
% hold on
% plot(Fv, abs(FXdcoc(Iv))*2)
% hold on
% plot(Fv, abs(FXdcoc2(Iv))*2)
% hold off
% grid
% title('Fourier Transform Of D-C Offset Corrected Signal ‘X’')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% [FXn_max,Iv_max] = max(abs(FXdcoc(Iv))*2);          % Get Maximum Amplitude, & Frequency Index Vector
% 
% % butterworth filter
% Wp = 2*Fv(Iv_max)/Fn;                               % Passband Frequency (Normalised)
% Ws = Wp*2;                                          % Stopband Frequency (Normalised)
% Rp = 10;                                            % Passband Ripple (dB)
% Rs = 30;                                            % Stopband Ripple (dB)
%[n,Wn] = buttord(Wp,Ws,Rp,Rs);                      % Butterworth Filter Order
[b,a] = butter(8,0.01);                               % Butterworth Transfer Function Coefficients
[b1,a1] = butter(8,0.001);                               % Butterworth Transfer Function Coefficients
% butterworth filter


% % bessel filter
% n1 = 8;                                               % number of poles (order)
% cut_off_freq = 20000;                                  % cut off frequency in MHz
% W0 = cut_off_freq*2*pi;
% [b,a] = besself(8,0.01);                               % bessel filter transfer function coefficients
% % bessel filter
r1 = 200;
rp = 35000;
cp = 1.8E-10;
cp1 = 2.2E-10;
cp2 = 2.4E-10;
cp3 = 2.6E-10;
cp4 = 2.8E-10;
V1 = 6.2*((r1/(r1+rp))+(rp/(r1+rp))*exp(-(t)*1E-6*(rp+r1)/(rp*r1*cp)));

[SOS,G] = tf2sos(b,a);                              % Convert to Second-Order-Section For Stability
% figure(3)
% freqz(SOS, 4096, Fs);                               % Filter Bode Plot
% title('Lowpass Filter Bode Plot')
S = filtfilt(SOS,G,X);
S2 = filtfilt(SOS,G,X2);                              % Filter ‘X’ To Recover ‘S’
figure(4)
%plot(t, X, 'LineWidth',0.5)                                          % Plot ‘X’
hold on
plot(t, X, '-r', 'LineWidth',0.5)                   % Plot ‘S’
hold on
plot(t, X2, '-g', 'LineWidth',0.5)                   % Plot ‘S’
hold on
plot(t, V1, '-b', 'LineWidth',0.5)                   % Plot ‘S’
hold off
ylim([-0.05 14])
grid
legend('Original Signal', 'Filtered Signal', 'Location','N')
title('Bubble signals')
xlabel('Time (\mus)')
ylabel('Voltage V1 (V)')

% fig = figure;
% left_color = [1 .5 0];
% right_color = [0 .5 0];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);

avgpt1 = 9;
avgpt2 = 16;

V = mean(X2((avgpt1*1E-6/(T)):(avgpt2*1E-6/(T))));
V1mean = mean(X((avgpt1*1E-6/(T)):(avgpt2*1E-6/(T))));
Vp = V-V1mean;



% hold on
% yyaxis left
% plot(t, X2, 'LineWidth',2.5)
% ylabel('Voltage V (V)')
% yyaxis right
% hold on
% plot(t, X, 'LineWidth',2.5)
% hold on
% plot(t, V1, '--b', 'LineWidth',2.5)
% hold off
% ylim([0 8])
% ylabel('Voltage V1 (V)')
% hold off
% set(gca,'FontSize',14)
% title('Oscilloscope signals for shunt voltage (Rs = 500 ohm)')
% xlabel('Time (\mus)')
% legend('Pulse Generator Voltage', 'Shunt Voltage', 'Fitting shunt voltage (Cp = 75 pF)','Location','N')


decimatingfactor = 1;
S1=decimate(S,decimatingfactor);
S3=decimate(S2,decimatingfactor);
t1 = decimate(t,decimatingfactor);




Iq1 = X/r1;
Iq = filtfilt(SOS,G,Iq1);
% Iq = Iq1;


[b2,a2] = butter(8,0.1); 
[SOS2,G2] = tf2sos(b2,a2); 

Iq2 = filtfilt(SOS2,G2,Iq1);



t1 =t1';
Iq =Iq*1e6;
Iq2 =Iq2*1e6;
% plot(t1,Iq)
% hold on
plot(t1,Iq2,'-r')
t1 = t1-104.78;


dt = t1(2)-t1(1);
delete('M1.csv');
delete('M2.csv');
delete('M3.csv');
delete('M4.csv');
delete('M5.csv');
delete('M6.csv');
delete('M7.csv');
delete('M8.csv');
delete('M9.csv');
delete('M10.csv');
delete('M11.csv');
delete('M12.csv');
delete('M13.csv');
delete('M14.csv');
delete('M15.csv');
delete('M16.csv');
delete('M17.csv');
delete('M18.csv');
delete('M19.csv');
delete('M20.csv');


M1=csvWrite(1,100,t1,Iq,Iq2);
M2=csvWrite(2,100,t1,Iq,Iq2);
M3=csvWrite(3,100,t1,Iq,Iq2);
M4=csvWrite(4,100,t1,Iq,Iq2);
M5=csvWrite(5,100,t1,Iq,Iq2);
M6=csvWrite(6,100,t1,Iq,Iq2);
M7=csvWrite(7,100,t1,Iq,Iq2);
M8=csvWrite(8,100,t1,Iq,Iq2);
M9=csvWrite(9,100,t1,Iq,Iq2);
M10=csvWrite(10,100,t1,Iq,Iq2);
M11=csvWrite(11,100,t1,Iq,Iq2);
M12=csvWrite(12,100,t1,Iq,Iq2);
M13=csvWrite(13,100,t1,Iq,Iq2);
M14=csvWrite(14,100,t1,Iq,Iq2);
M15=csvWrite(15,100,t1,Iq,Iq2);
M16=csvWrite(16,100,t1,Iq,Iq2);
M17=csvWrite(17,100,t1,Iq,Iq2);
M18=csvWrite(18,100,t1,Iq,Iq2);
M19=csvWrite(19,100,t1,Iq,Iq2);
M20=csvWrite(20,100,t1,Iq,Iq2);

figure
plot(M2(:,1), M2(:,2), 'LineWidth',1.5)
hold on
p1 = plot(M2(:,1), M2(:,3),'-r');
p1.Color(4) = 0.25;
ylim([-500 3000])
hold off





% legend('Shunt Voltage', 'Pore Conductance (cp = 200pF)', 'Pore Conductance (cp = 220pF)', 'Pore Conductance (cp = 240pF)', 'Pore Conductance (cp = 260pF)', 'Pore Conductance (cp = 280pF)','Location','N')
% plot(t, S, '-r', 'LineWidth',1.5)
% grid
% title('Bubble signals')
% xlabel('Time (\mus)')
% ylabel('Voltage V1 (V)')
