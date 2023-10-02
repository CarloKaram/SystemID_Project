close all
clear all

load PWM_left
load PWM_right
load disp_x
load disp_y
load disp_theta
load Vel_x
load Vel_y
load Vel_ang

%% Data Preprocessing : handling nonzero means

u_1 = u_1 - (1/1000)*sum(u_1) ;
u_2 = u_2 - (1/1000)*sum(u_2) ;

LinearX_both = LinearX_both - (1/1000)*sum(LinearX_both) ;
Angluar_both = Angluar_both - (1/1000)*sum(Angluar_both) ;

figure(1)
autocorr(LinearX_both)
title('Autocorrelation of Input into Left Motor')

figure(2)
autocorr(Angluar_both)
title('Autocorrelation of Input into Right Motor')

%% Spectral Estimates of Inputs and Outputs

Nfft = 512;

figure(1)
Nwin = 32;
[Pxx, w] = pwelch(u_1,Nwin,Nwin/2,Nfft);
plot(w,Pxx); 
title('Spectral Estimate of PWM input into left motor')
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Amplitude')
xlim([0 3.142])
grid on

figure(2)
Nwin = 32;
[Pxx, w] = pwelch(u_2,Nwin,Nwin/2,Nfft);
plot(w,Pxx); 
title('Spectral Estimate of PWM input into right motor')
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Amplitude')
xlim([0 3.142])
grid on

figure(3)
Nwin = 32;
[Pxx, w] = pwelch(LinearX_both,Nwin,Nwin/2,Nfft);
plot(w,Pxx); 
title('Spectral Estimate of Linear Velocity')
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Amplitude')
xlim([0 3.142])
grid on

figure(4)
Nwin = 32;
[Pxx, w] = pwelch(Angluar_both,Nwin,Nwin/2,Nfft);
plot(w,Pxx); 
title('Spectral Estimate of Angular Velocity')
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Amplitude')
xlim([0 3.142])
grid on

%% Coherence Plots

figure(1)
[Cu_V, w] = mscohere(u_1,LinearX_both,50,25);
plot(w,Cu_V)
grid on
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Magnitude')
xlim([0 3.142])
title('Coherence Estimate between Linear Velocity and PWM input to left motor')

figure(2)
[Cu_omega, w] = mscohere(u_1,Angluar_both,50,25);
plot(w,Cu_omega)
grid on
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Magnitude')
xlim([0 3.142])
title('Coherence Estimate between Angular Velocity and PWM input to left motor')

figure(3)
[Cu_V, w] = mscohere(u_2,LinearX_both,50,25);
plot(w,Cu_V)
grid on
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Magnitude')
xlim([0 3.142])
title('Coherence Estimate between Linear Velocity and PWM input to right motor')

figure(4)
[Cu_omega, w] = mscohere(u_2,Angluar_both,128,64);
plot(w,Cu_omega)
grid on
xlabel('Normalized Frequencies (rad/sample)')
ylabel('Magnitude')
xlim([0 3.142])
title('Coherence Estimate between Angular Velocity and PWM input to right motor')

%% Cross Correlations

figure(1)
crosscorr(u_1,LinearX_both)
title('Cross-correlation plot between PWM input of left motor and Linear Velocity')

figure(2)
crosscorr(u_2,LinearX_both)
title('Cross-correlation plot between PWM input of right motor and Linear Velocity')

figure(3)
crosscorr(u_1,Angluar_both)
title('Cross-correlation plot between PWM input of left motor and Angular Velocity')

figure(4)
crosscorr(u_2,Angluar_both)
title('Cross-correlation plot between PWM input of right motor and Angular Velocity')
