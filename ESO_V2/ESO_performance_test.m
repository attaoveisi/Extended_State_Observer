clear all
close all
%% System model
load('identified_n_8.mat');
load('S.mat');
H = B(:,4);
% H2 = H1;
B = B(:,1:2);
D = D(:,1:2);
nx = size(A,1);
nu = size(B,2);
nw = size(H,2);
ny = size(C,1);
%% Disturbance model
nv = size(S,1);
nh = nv;
F = eye(nv);
E = zeros(nw,nv);
for i = 1:6
E(:,2*i-1:2*i) = [1 1];
end
%% ESO gains
load('La.mat')
L1 = La(1:nx,:);
L2 = La(nx+1:nx+nv,:);
%% Simulation
Tsim = 20;%Simulation time (s)
sim('ESO_test')
figure(1)
subplot(3,1,1)
plot(t,w)
ylabel('$w$','interpreter','latex','fontsize',14)
subplot(3,1,2)
plot(t,w_hat,'r')
ylabel('$\hat{w}$','interpreter','latex','fontsize',14)
subplot(3,1,3)
plot(t,w-w_hat,'g')
xlabel('Time (s)','interpreter','latex','fontsize',14)
ylabel('$w-\hat{w}$','interpreter','latex','fontsize',14)
% ss_index = t>=8 & t<=t(end);
% axes('Position',[0.7 0.23 0.2 0.2]);
% plot(t(ss_index),w(ss_index)-w_hat(ss_index),'g')
figure(2)
subplot(4,1,1)
plot(t,x(:,1),t,x_hat(:,1),'r')
hold on
plot(t,x(:,1)-x_hat(:,1),'g')
legend('True','Estimated','Error')
subplot(4,1,2)
plot(t,x(:,2),t,x_hat(:,2),'r')
hold on
plot(t,x(:,2)-x_hat(:,2),'g')
subplot(4,1,3)
plot(t,x(:,3),t,x_hat(:,3),'r')
hold on
plot(t,x(:,3)-x_hat(:,3),'g')
subplot(4,1,4)
plot(t,x(:,4),t,x_hat(:,4),'r')
hold on
plot(t,x(:,4)-x_hat(:,4),'g')
%
figure(3)
subplot(4,1,1)
plot(t,x(:,5),t,x_hat(:,5),'r')
hold on
plot(t,x(:,5)-x_hat(:,5),'g')
legend('True','Estimated','Error')
subplot(4,1,2)
plot(t,x(:,6),t,x_hat(:,6),'r')
hold on
plot(t,x(:,6)-x_hat(:,6),'g')
subplot(4,1,3)
plot(t,x(:,7),t,x_hat(:,7),'r')
hold on
plot(t,x(:,7)-x_hat(:,7),'g')
subplot(4,1,4)
plot(t,x(:,8),t,x_hat(:,8),'r')
hold on
plot(t,x(:,8)-x_hat(:,8),'g')