clear all
%% 
load('identified_modal_250.mat');
load('S.mat');
S=S(1:2,1:2);
n=size(A,1);
n_w=size(S,1);
n_y=size(C,1);
%
Cw=zeros(size(H,2),n_w);
Cw=[1 0];
H=H*Cw;
Q=eye(size(H));
n_xi=size(Q,2);
n_t=n+n_w+n_xi;
A0=zeros(n_t,n_t);
B0=zeros(n_t,n_xi);
C0=zeros(n_y,n_t);
%
A0(1:n,:)=[A H Q];
A0(n+1:n+n_w,n+1:n+n_w)=S;
B0(n+n_w+1:n_t,:)=eye(n_xi);
C0(:,1:n)=C;
Cz=C0;
%% Sector boundeds of the gain function [a,b]
a=0.5;
b=1;
%% LMI optimization
c1=1e+09;
c2=1e+07;
threshold=1e-06;%Threshold for stric inequalities
cvx_begin sdp
variable  P(n_t,n_t) symmetric
variables Y(n_t,n_y) gam(1,1) tau(1,1)
%
psi11=P*A0+A0'*P+(Cz'*Cz)-tau*a*b*(C0'*C0);
psi12=-Y+(tau/2)*(a+b)*C0';
psi13=P*B0;
psi22=-tau*eye(n_y);
psi23=zeros(n_y,n_xi);
psi33=-gam*eye(n_xi);
Psi=[psi11  psi12  psi13;
     psi12' psi22  psi23;
     psi13' psi23' psi33];
%
minimize(gam)
subject to
gam>=threshold;
tau>=threshold;
% P-c1^(-1)*eye(n_t)>=0;
% [-c2^2*eye(n_y) Y';
%   Y  -eye(n_t)]<=0;
Psi<=-threshold*eye(size(Psi));
cvx_end
%% Results
gam_opt=gam^0.5;
L=P\Y;