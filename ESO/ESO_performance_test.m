clear all
close all
%% System model
% load('identified_modal.mat');
% Am = A;
% Bm = B;
% Cm = C;
% Dm = D(:,1:3);
% Hm = H;
% clear A B C D H
%% ESO model
load('identified_modal_250.mat');
D = D(:,1:3);
Am = A;
Bm = B;
Cm = C;
Dm = D;
Hm = H;
load('S.mat');
S = blkdiag(S,zeros(1,1));
n_x = size(A,1);
n_xw = size(S,1);
n_w = size(H,2);
n_h = 1;
E = zeros(n_xw,n_h);
E(n_xw,:)=1;
Cw = zeros(n_w,n_xw);
for i=1:floor(n_xw/2)
Cw(:,2*i-1:2*i) = [1 1];
end
Cw(:,n_xw) = 1;
load('ESO.mat')
L1 = L(1:n_x,:);
L2 = L(n_x+1:n_x+n_xw,:);
%% Simulation
Tsim = 10;%Simulation time (s)
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
ss_index = t>=8 & t<=t(end);
axes('Position',[0.7 0.23 0.2 0.2]);
plot(t(ss_index),w(ss_index)-w_hat(ss_index),'g')