%% EN2570 - FIR Bandstop Filter Design
% By K.H.M.K.D.Kehelella - 170302T 

clc;
clear all;
close all;
%% Required Filter Specifications

A = 3;
B = 0;
C = 2;

A_p = 0.03+(0.01*A); % dB %%max passband ripple
A_a = 45+B; %dB %%min stopband attenuation
Omega_p1 = (C*100)+400; %rad/s %%lower passband edge
Omega_p2 = (C*100)+950; %rad/s %%upper passband edge
Omega_a1 = (C*100)+500; %rad/s %%lower stopband edge
Omega_a2 = (C*100)+800; %rad/s %%upper stopband edge
Omega_s = 2*(C*100+1300); %rad/s %%sampling freqency

%% Derived Filter Specifications

Bt1 = Omega_a1-Omega_p1; %rad/s %%lower transition width
Bt2 = Omega_p2-Omega_a2;  %rad/s %%upper transisiton width
Bt = min(Bt1,Bt2); %rad/s %%critical transition width
Omega_c1 = Omega_p1+Bt/2; %rad/s %%lower cutoff frequency
Omega_c2 = Omega_p2-Bt/2; %rad/s %%upper cutoff frequency
T = 2*pi/Omega_s; %s %%sampling period

%% Kaiser Window Parameters

% calculating delta
d_P = (10^(0.05*A_p) - 1)/ (10^(0.05*A_p) + 1); 
d_A = 10^(-0.05*A_a);
delta = min(d_P,d_A);

Aa = -20*log10(delta);  % Actual stopband attenuation
Ap = 20*log10(1+delta/1-delta); % Actual passband ripple

% Calculating alpha
if Aa<=21                
    alpha = 0;
elseif Aa>21 && Aa<= 50
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa-8.7);
end

% Calculating D
if Aa <= 21             
    D = 0.9222;
else
    D = (Aa-7.95)/14.36;
end

% Calculating order of the filter N
N = ceil(Omega_s*D/Bt +1); 
if mod(N,2) == 0 %If lowest N is even then this will give lowest odd value otherwise it s already lowest odd value
    N = N+1;
end

% Length of the filter
n = -(N-1)/2:1:(N-1)/2;   

% Calculating beta
beta = alpha*sqrt(1-(2*n/(N-1)).^2);

%% Generating Io_(alpha)

bessellimit = 100;

Io_alpha = 1;
for k = 1:bessellimit
    val_k = (1/factorial(k)*(alpha/2).^k).^2;
    Io_alpha = Io_alpha + val_k;
end

%% Generating Io_(beta) 

Io_beta = 1;
for m = 1:bessellimit
    val_m = (1/factorial(m)*(beta/2).^m).^2;
    Io_beta = Io_beta + val_m;
end

%% Obtaining Kaiser Window w_k(nT)

wk_nT = Io_beta/Io_alpha;

figure
stem(n,wk_nT)
xlabel('n')
ylabel('Amplitude')
title('Kaiser Window in Time Domain');

%% Generating Impulse Response h(nT)

n_L= -(N-1)/2:1:-1;
hnt_L = 1./(n_L*pi).*(sin(Omega_c1*n_L*T)-sin(Omega_c2*n_L*T));

n_R = 1:1:(N-1)/2;
hnt_R = 1./(n_R*pi).*(sin(Omega_c1*n_R*T)-sin(Omega_c2*n_R*T));

hnt_0 = 1+ 2/Omega_s*(Omega_c1-Omega_c2);

n = [n_L,0,n_R];
h_nT = [hnt_L,hnt_0,hnt_R];

%% Applying the Kaiser Window to the filter 

Hw_nT = h_nT.*wk_nT;

%% Plotting the causal impulse response

n_shifted = [0:1:N-1]; %we must shift this much of samples to obtain causal impulse response
figure
stem(n_shifted,Hw_nT); axis tight;
xlabel('n')
ylabel('Amplitude')
title(strcat(['Kaiser Window Filter Causal Impulse Response in Time Domain']));

%% Plotting the magnitude response of filter in the range (0,Os/2)

fvtool(Hw_nT);
figure
[Hw,f] = freqz(Hw_nT);
w_1 = f*Omega_s/(2*pi);
log_Hw = 20*log10(abs(Hw));
plot(w_1,log_Hw)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat(['Kaiser Window Filter Magnitude Response in Frequency Domain']));

%% Plotting the magnitude response of filter in Lower Passband 

figure
finish = round((length(w_1)/(Omega_s/2)*Omega_c1));
wpass_l = w_1(1:finish);
hpass_l = log_Hw(1:finish);
plot(wpass_l,hpass_l)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Filter Lower Passband in Frequency Domain');

%% Plotting the magnitude response of filter in Upper Passband 

figure
start = round(length(w_1)/(Omega_s/2)*Omega_c2);
wpass_h = w_1(start:length(w_1));
hpass_h = log_Hw(start:length(w_1));
plot(wpass_h,hpass_h)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Filter Upper Passband in Frequency Domain');

%% Plotting the Rectangular Window filter response

figure
stem(n_shifted,h_nT); axis tight;
xlabel('n')
ylabel('Amplitude')
title(strcat(['Rectangular window Filter Response in Time Domain']));

figure
[hw,f] = freqz(h_nT);
w_2 = f*Omega_s/(2*pi);
log_H = 20*log10(hw);
plot(w_2,log_H)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat(['Rectangular Window Filter Response in Frequency Domain']));

%% Input signal generation 

% component frequencies of the input
Omega_1 = Omega_c1/2;                    
Omega_2 = Omega_c1 + (Omega_c2-Omega_c1)/2;
Omega_3 = Omega_c2 + (Omega_s/2-Omega_c2)/2;

% generate discrete signal and evelope
samples = 500;
n1 = 0:1:samples;                
                      
X_nT = cos(Omega_1.*n1.*T)+cos(Omega_2.*n1.*T)+cos(Omega_3.*n1.*T); 

%% Using DFT to check the filtering 

% Filtering using frequency domain multiplication 

len_fft = length(X_nT)+length(Hw_nT)-1; % length for fft in x dimension
x_fft = fft(X_nT,len_fft);
Hw_nT_fft = fft(Hw_nT,len_fft);
out_fft = Hw_nT_fft.*x_fft; % A shift in time is added here
out = ifft(out_fft,len_fft);
out = out(1:length(out)-floor(N)+1);
rec_out = out(floor(N/2)+1:length(out)-floor(N/2));% account for shifting delay

% Ideal Output Signal
ideal_out = cos(Omega_1.*n1.*T)+cos(Omega_3.*n1.*T);

%% Plots of final results

% Time domain representation of input signal before filtering

figure 
subplot(2,1,1)
stem(n1,X_nT)
xlabel('n')
ylabel('Amplitude')
title(strcat(['Input signal',' ','in Time Domain']));

% Time domain representation of output signal after filtering

subplot(2,1,2)
stem(n1,out)
xlabel('n')
ylabel('Amplitude')
title(strcat(['Output signal',' ','in Time Domain']));

% Frequency domain representation of input signal before filtering

figure 
subplot(2,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
x_fft = fft(X_nT,len_fft);
x_fft_plot = [abs([x_fft(len_fft/2+1:len_fft)]),abs(x_fft(1)),abs(x_fft(2:len_fft/2+1))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_plot); axis tight;
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Input signal',' ','in Frequency Domain']));

% Frequency domain representation of output signal after filtering

subplot(2,1,2)
len_fft = 2^nextpow2(numel(n1))-1;
xfft_out = fft(rec_out,len_fft);
x_fft_out_plot = [abs([xfft_out(len_fft/2+1:len_fft)]),abs(xfft_out(1)),abs(xfft_out(2:len_fft/2+1))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_out_plot); axis tight;
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Output signal',' ','in Frequency Domain']));

% Comparing Frequency domain representation of ideal output signal with frequency domain representation of output signal after filtering

figure 
subplot(2,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
x_fft = fft(ideal_out,len_fft);
x_fft_ideal_out_plot = [abs([x_fft(len_fft/2+1:len_fft)]),abs(x_fft(1)),abs(x_fft(2:len_fft/2+1))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_ideal_out_plot); axis tight;
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Ideal output signal',' ','in Frequency Domain']));

% Frequency domain representation of output signal after filtering to compair
subplot(2,1,2)
len_fft = 2^nextpow2(numel(n1))-1;
xfft_out = fft(rec_out,len_fft);
x_fft_out_plot = [abs([xfft_out(len_fft/2+1:len_fft)]),abs(xfft_out(1)),abs(xfft_out(2:len_fft/2+1))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_out_plot); axis tight;
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Output signal',' ','in Frequency Domain']));