clear all
clc
%% Measured Values

load Vel_x ;
load Vel_ang ;
load PWM_left
load PWM_right


eLin = sqrt((var(LinearX_both))/5)*randn(1,1000) ;
eAng = sqrt((var(Angluar_both)/5))*randn(1,1000) ;

LinearVel = LinearX_both - mean(LinearX_both) + eLin ;
AngularVel = Angluar_both - mean(Angluar_both) + eAng ;

%% Input Generation



%% Data Generation
for k = 3:1000
    y1(k) = (-0.3395*u_1(k-2)+0.3681*u_2(k-2))  ;
    y2(k) = (-0.0699*u_1(k-1)+0.0758*u_2(k-1))  ;
    %y1(k) = (-0.3354*u_1(k-2)+0.3636*u_2(k-2))  ;
    %y2(k) = (-0.0656*u_1(k-2)+0.0711*u_2(k-2))  ;
end

resLin = abs(LinearVel(1:end) - y1) ;
resAng = abs(AngularVel(1:end) - y2) ;

t = 0:0.1:100 ;

figure(1)
plot(t(1:end-1),LinearVel(1:end),'-r')
grid on;
hold on;
plot(t(1:end-1),y1,'--b')
xlim([0 10])
title('Actual vs. Estimated Linear Velocity (noisy)')
legend('Actual Velocity','Estimated Velocity')

figure(2)
plot(t(1:end-1),AngularVel(1:end),'-r')
grid on;
hold on;
plot(t(1:end-1),y2,'--b')
xlim([0 10])
title('Actual vs. Estimated Angular Velocity (noisy)')
legend('Actual Velocity','Estimated Velocity')

figure(3)
autocorr(resLin)

figure(4)
autocorr(resAng)

figure(5)
crosscorr(u_1(1:end-2),resLin)
xlim([-20 20])
ylim([-1 1])

figure(6)
crosscorr(u_2(1:end-2),resLin) ;
xlim([-20 20])
ylim([-1 1])

figure(7)
crosscorr(u_1(1:end-2),resAng) ;
xlim([-20 20])
ylim([-1 1])

figure(8)
crosscorr(u_2(1:end-2),resAng) ;
xlim([-20 20])
ylim([-1 1])
