
close all
clear all                                                               
clc   

%% initialization

% initial conditions                                                    
x0 = [5;0]; 
xHat0 = [10;2];
          
% physical variables 
gamma = -9.81;  % gravity constant
lambda = 0.8;   % restitution coefficent
rho = 0.1; % friction coefficient

% simulation horizon                                                    
T = 15;                                                                 
J = 20;                                                                 
                                                                        
% rule for jumps                                                        
% rule = 1 -> priority for jumps                                        
% rule = 2 -> priority for flows                                        
% rule = 3 -> no priority, random selection when simultaneous conditions
rule = 1;                                                               
                                                                        
%solver tolerances
RelTol = 1e-6;
MaxStep = 1e-3;

% Observer 
Ac = [0,1;0,-rho];
Ad = [-1,0;0,-lambda];
C = [1,0];

% gain for update at jumps only
Lc =0;
Ld = [-1;-0.1085]; %yalmip lambda = 0.8 rho = 0.1

% eigenvaues of the sampled system (to check that <1 on [tau_m,tau_M]
tau = linspace(0,10,1000);
max_eig = zeros(length(tau),2);
for ind=1:length(tau)
    X = eig(expm((Ac-Lc*C)*tau(ind))*(Ad-Ld*C));
    max_eig(ind,:)=abs(X)';
end
figure(1)
plot(tau,max_eig)
title('Eigenvalues of $\exp((A_c-L_cH_c)\tau)(A_d-L_dH_d)$','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')


%% Simulation

sim('BouncingBall');

%% Plots

% construction of jump vector common to x and xHat to plot error
jRes = zeros(size(j));
for ind=2:length(jRes)
    if j(ind)~=j(ind-1) || jHat(ind)~=jHat(ind-1)
        jRes(ind) = jRes(ind-1)+1;
    else 
        jRes(ind) = jRes(ind-1);
    end
end

modificatorF{1} = '-';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 2;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.2;

% plot solution
figure(2) 
clf
subplot(2,1,1), plotHarc([t,tHat],[j,jHat],[x(:,1),xHat(:,1)],[],modificatorF,modificatorJ);
grid on
leg1 = legend('$\theta$','$\hat{\theta}$');
set(leg1,'Interpreter','latex','Fontsize',20)
xlabel('$t$ [s]','Fontsize',15,'Interpreter','latex')
subplot(2,1,2), plotHarc([t,tHat],[j,jHat],[x(:,2),xHat(:,2)],[],modificatorF,modificatorJ);
grid on
leg2 = legend('$\omega$','$\hat{\omega}$');
set(leg2,'Interpreter','latex','Fontsize',20)
xlabel('$t$ [s]','Fontsize',15,'Interpreter','latex')

% plot error
figure(3) 
clf
subplot(1,1,1)
plotHarc([t,t],[jRes,jRes],[-x(:,1)+xHat(:,1),-x(:,2)+xHat(:,2)],[],modificatorF,modificatorJ)
leg=legend('$\hat{\theta}-\theta$','$\hat{\omega}-\omega$');
set(leg, 'Interpreter', 'latex','Fontsize',20)
xlabel('$t$ [s]','Fontsize',15,'Interpreter','latex')
grid on

