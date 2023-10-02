%% Load variables

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

%% Subtract mean to remove offset and transpose to get column vectors

u_1 = u_1-mean(u_1);
u_2 = u_2-mean(u_2);
x_both = (x_both-mean(x_both))';
y_both = (y_both-mean(y_both))';
theta_both = (theta_both-mean(theta_both))';
LinearX_both = (LinearX_both-mean(LinearX_both))';
LinearY_both = (LinearY_both-mean(LinearY_both))';
Angluar_both = (Angluar_both-mean(Angluar_both))';

% Set noise (0 or 1) and SNR 
noise = 1;
SNR = 5;

N = length(u_1);

%% Vx (ARX - increasing order models)

var_Vx = var(LinearX_both);
e_Vx = noise*sqrt(var_Vx/SNR)*randn(length(LinearX_both),1);
% SNR = power_signal/power_noise = 5

LinearX_both = LinearX_both + e_Vx;

THETA_Vx = NaN(12,4);
PHI_Vx = [];
EPS_Vx = [];
V_Vx = zeros(1,4);
sigma_e = zeros(1,4);
COV = [];

for ord=2:5
    phi0 = [];
    
    for k = ord-1:-1:1
        phi0 = [phi0 -LinearX_both(k:end+k-ord) u_1(k:end+k-ord) u_2(k:end+k-ord)];
    end
    
    theta0 = inv(phi0'*phi0)*phi0'*LinearX_both(ord:end);
    
    THETA_Vx(1:length(theta0),ord-1) = theta0;
    phi_aug = [phi0;zeros(ord-2,3*(ord-1))];
    PHI_Vx =[PHI_Vx phi_aug]; 
    
    eps0 = LinearX_both(ord:end) - phi0*theta0;
    eps_aug = [eps0;zeros((ord-2),1)];
    EPS_Vx = [EPS_Vx eps_aug];
    
    V_Vx(ord-1) = 0.5*EPS_Vx(:,ord-1)'*EPS_Vx(:,ord-1);
    sigma_e(ord-1) = 2*V_Vx(ord-1)/(N-3*(ord-1));
    
    cov0 = sigma_e(ord-1)*inv(phi0'*phi0);
    cov_aug = [cov0;zeros(9-3*(ord-2),3*(ord-1))];
    COV = [COV cov_aug];
    
end

figure(1)
autocorr(EPS_Vx(:,1))

figure(2)
autocorr(EPS_Vx(1:end-1,2))

figure(3)
autocorr(EPS_Vx(1:end-2,3))

figure(4)
autocorr(EPS_Vx(1:end-3,4))


figure(5)
crosscorr(PHI_Vx(:,1),EPS_Vx(:,1))

figure(6)
crosscorr(PHI_Vx(:,2),EPS_Vx(:,1))

figure(7)
crosscorr(PHI_Vx(:,3),EPS_Vx(:,1))


figure(8)
crosscorr(PHI_Vx(1:end-1,4),EPS_Vx(1:end-1,2))

figure(9)
crosscorr(PHI_Vx(1:end-1,5),EPS_Vx(1:end-1,2))

figure(10)
crosscorr(PHI_Vx(1:end-1,6),EPS_Vx(1:end-1,2))

figure(11)
crosscorr(PHI_Vx(1:end-1,7),EPS_Vx(1:end-1,2))

figure(12)
crosscorr(PHI_Vx(1:end-1,8),EPS_Vx(1:end-1,2))

figure(13)
crosscorr(PHI_Vx(1:end-1,9),EPS_Vx(1:end-1,2))


%% w (ARX - increasing order models)

var_w = var(Angluar_both);
e_w = noise*sqrt(var_w/SNR)*randn(length(Angluar_both),1);
% SNR = power_signal/power_noise = 5

Angluar_both = Angluar_both + e_w;

THETA_w = NaN(12,4);
PHI_w = [];
EPS_w = [];
V_w = zeros(1,4);
sigma_e = zeros(1,4);
COV = [];

for ord=2:5
    phi0 = [];
    
    for k = ord-1:-1:1
        phi0 = [phi0 -Angluar_both(k:end+k-ord) u_1(k:end+k-ord) u_2(k:end+k-ord)];
    end
    
    theta0 = inv(phi0'*phi0)*phi0'*Angluar_both(ord:end);
    
    THETA_w(1:length(theta0),ord-1) = theta0;
    phi_aug = [phi0;zeros(ord-2,3*(ord-1))];
    PHI_w =[PHI_w phi_aug]; 
    
    eps0 = Angluar_both(ord:end) - phi0*theta0;
    eps_aug = [eps0;zeros((ord-2),1)];
    EPS_w = [EPS_w eps_aug];
    
    V_w(ord-1) = 0.5*EPS_w(:,ord-1)'*EPS_w(:,ord-1);
    sigma_e(ord-1) = 2*V_w(ord-1)/(N-3*(ord-1));
    
    cov0 = sigma_e(ord-1)*inv(phi0'*phi0);
    cov_aug = [cov0;zeros(9-3*(ord-2),3*(ord-1))];
    COV = [COV cov_aug];
    
end

figure(1)
autocorr(EPS_w(:,1))

figure(2)
autocorr(EPS_w(1:end-1,2))

figure(3)
autocorr(EPS_w(1:end-2,3))

figure(4)
autocorr(EPS_w(1:end-3,4))


figure(5)
crosscorr(PHI_w(:,1),EPS_w(:,1))

figure(6)
crosscorr(PHI_w(:,2),EPS_w(:,1))

figure(7)
crosscorr(PHI_w(:,3),EPS_w(:,1))


figure(8)
crosscorr(PHI_w(1:end-1,4),EPS_w(1:end-1,2))

figure(9)
crosscorr(PHI_w(1:end-1,5),EPS_w(1:end-1,2))

figure(10)
crosscorr(PHI_w(1:end-1,6),EPS_w(1:end-1,2))

figure(11)
crosscorr(PHI_w(1:end-1,7),EPS_w(1:end-1,2))

figure(12)
crosscorr(PHI_w(1:end-1,8),EPS_w(1:end-1,2))

figure(13)
crosscorr(PHI_w(1:end-1,9),EPS_w(1:end-1,2))


%% Vx (ideal FIR - inputs 2 steps before)

var_Vx = var(LinearX_both);
e_Vx = noise*sqrt(var_Vx/SNR)*randn(length(LinearX_both),1);
% SNR = power_signal/power_noise = 5

LinearX_both = LinearX_both + e_Vx;

figure(1)
crosscorr(u_1,LinearX_both)
figure(2)
crosscorr(u_2,LinearX_both)
figure(3)
autocorr(LinearX_both)

phi = [u_1(1:end-2) u_2(1:end-2)];
y = LinearX_both(3:end);
theta = inv(phi'*phi)*phi'*y;
eps = y - phi*theta;

V = 0.5*eps'*eps;
sigma = 2*V/(N-2);
cov = sigma*inv(phi'*phi);

figure(4)
autocorr(eps)
figure(5)
crosscorr(phi(:,1),eps)
figure(6)
crosscorr(phi(:,2),eps)

%% w (ideal FIR - inputs 2 steps before)

var_w = var(Angluar_both);
e_w = noise*sqrt(var_w/SNR)*randn(length(Angluar_both),1);
% SNR = power_signal/power_noise = 5

Angluar_both = Angluar_both + e_w;

figure(1)
crosscorr(u_1,Angluar_both)
figure(2)
crosscorr(u_2,Angluar_both)
figure(3)
autocorr(Angluar_both)

phi = [u_1(1:end-2) u_2(1:end-2)];
y = Angluar_both(3:end);
theta = inv(phi'*phi)*phi'*y;
eps = y - phi*theta;

V = 0.5*eps'*eps;
sigma = 2*V/(N-2);
cov = sigma*inv(phi'*phi);

figure(4)
autocorr(eps)
figure(5)
crosscorr(phi(:,1),eps)
figure(6)
crosscorr(phi(:,2),eps)


%% Vx (FIR - increasing order models)

var_Vx = var(LinearX_both);
e_Vx = noise*sqrt(var_Vx/SNR)*randn(length(LinearX_both),1);
% SNR = power_signal/power_noise = 5

LinearX_both = LinearX_both + e_Vx;

THETA_Vx = NaN(8,4);
PHI_Vx = [];
EPS_Vx = [];
V_Vx = zeros(1,4);
sigma_e = zeros(1,4);
COV = [];

for ord=2:5
    phi0 = [];
    
    for k = ord-1:-1:1
        phi0 = [phi0 u_1(k:end+k-ord) u_2(k:end+k-ord)];
    end
    
    theta0 = inv(phi0'*phi0)*phi0'*LinearX_both(ord:end);
    
    THETA_Vx(1:length(theta0),ord-1) = theta0;
    phi_aug = [phi0;zeros(ord-2,2*(ord-1))];
    PHI_Vx =[PHI_Vx phi_aug]; 
    
    eps0 = LinearX_both(ord:end) - phi0*theta0;
    eps_aug = [eps0;zeros((ord-2),1)];
    EPS_Vx = [EPS_Vx eps_aug];
    
    V_Vx(ord-1) = 0.5*EPS_Vx(:,ord-1)'*EPS_Vx(:,ord-1);
    sigma_e(ord-1) = 2*V_Vx(ord-1)/(N-2*(ord-1));
    
    cov0 = sigma_e(ord-1)*inv(phi0'*phi0);
    cov_aug = [cov0;zeros(6-2*(ord-2),2*(ord-1))];
    COV = [COV cov_aug];
    
end

figure(1)
autocorr(EPS_Vx(:,1))

figure(2)
autocorr(EPS_Vx(1:end-1,2))

figure(3)
autocorr(EPS_Vx(1:end-2,3))

figure(4)
autocorr(EPS_Vx(1:end-3,4))


figure(5)
crosscorr(PHI_Vx(:,1),EPS_Vx(:,1))

figure(6)
crosscorr(PHI_Vx(:,2),EPS_Vx(:,1))


figure(7)
crosscorr(PHI_Vx(1:end-1,3),EPS_Vx(1:end-1,2))

figure(8)
crosscorr(PHI_Vx(1:end-1,4),EPS_Vx(1:end-1,2))

figure(9)
crosscorr(PHI_Vx(1:end-1,5),EPS_Vx(1:end-1,2))

figure(10)
crosscorr(PHI_Vx(1:end-1,6),EPS_Vx(1:end-1,2))


%% w (FIR - increasing order models)

var_w = var(Angluar_both);
e_w = noise*sqrt(var_w/SNR)*randn(length(Angluar_both),1);
% SNR = power_signal/power_noise = 5

Angluar_both = Angluar_both + e_w;

THETA_w = NaN(8,4);
PHI_w = [];
EPS_w = [];
V_w = zeros(1,4);
sigma_e = zeros(1,4);
COV = [];

for ord=2:5
    phi0 = [];
    
    for k = ord-1:-1:1
        phi0 = [phi0 u_1(k:end+k-ord) u_2(k:end+k-ord)];
    end
    
    theta0 = inv(phi0'*phi0)*phi0'*Angluar_both(ord:end);
    
    THETA_w(1:length(theta0),ord-1) = theta0;
    phi_aug = [phi0;zeros(ord-2,2*(ord-1))];
    PHI_w =[PHI_w phi_aug]; 
    
    eps0 = Angluar_both(ord:end) - phi0*theta0;
    eps_aug = [eps0;zeros((ord-2),1)];
    EPS_w = [EPS_w eps_aug];
    
    V_w(ord-1) = 0.5*EPS_w(:,ord-1)'*EPS_w(:,ord-1);
    sigma_e(ord-1) = 2*V_w(ord-1)/(N-2*(ord-1));
    
    cov0 = sigma_e(ord-1)*inv(phi0'*phi0);
    cov_aug = [cov0;zeros(6-2*(ord-2),2*(ord-1))];
    COV = [COV cov_aug];
    
end

figure(1)
autocorr(EPS_w(:,1))

figure(2)
autocorr(EPS_w(1:end-1,2))

figure(3)
autocorr(EPS_w(1:end-2,3))

figure(4)
autocorr(EPS_w(1:end-3,4))


figure(5)
crosscorr(PHI_w(:,1),EPS_w(:,1))

figure(6)
crosscorr(PHI_w(:,2),EPS_w(:,1))


figure(7)
crosscorr(PHI_w(1:end-1,3),EPS_w(1:end-1,2))

figure(8)
crosscorr(PHI_w(1:end-1,4),EPS_w(1:end-1,2))

figure(9)
crosscorr(PHI_w(1:end-1,5),EPS_w(1:end-1,2))

figure(10)
crosscorr(PHI_w(1:end-1,6),EPS_w(1:end-1,2))

