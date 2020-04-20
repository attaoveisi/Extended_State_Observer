clear all
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
%% Checking the observability
PBH(A',C')
PBH(S',E')
%% Augmented error model
Aa = [A H*E;
      zeros(nv,nx) S];
Ba = [zeros(nx,nh);F];
Ca = [C zeros(ny,nv)];
Cd = Ca; %Desired ouput
nxa = nx+nv;
%% Sector condition
alpha = 0.8;
beta = 1;
%% LMIs
%Inputs
kappa0 = 1e-01;
c1 = 1e+03;
c2 = 1e+04;
%The optimal value of muu without gain constraint on La is 0.0537
%The corresponding La norm is 2.8575e+06
% muu = 1;
cvx_begin sdp
variable P(nxa,nxa) symmetric
variable Y(nxa,ny)
variables tau(1,1) muu(1,1)
minimize(muu)
subject to
tau>=1e-08;
muu>=1e-08;
Pi_11 = P*Aa+Aa'*P+kappa0*eye(nxa)+(Cd'*Cd)-tau*alpha*beta*(Ca'*Ca);
Pi_12 = -Y+0.5*(alpha+beta)*tau*Ca';
Pi_13 = P*Ba;
Pi_22 = -tau*eye(ny);
Pi_23 = zeros(ny,nh);
Pi_33 = -muu*eye(nh);
[Pi_11 Pi_12 Pi_13;
 Pi_12' Pi_22 Pi_23;
 Pi_13' Pi_23' Pi_33]<=0;
P-c1^(-1)*eye(nxa)>=0;
[-c2^2*eye(ny) Y';
    Y -eye(nxa)]<=0;
cvx_end
La = P\Y;
save('La.mat');