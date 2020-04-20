clear all
%% Loading the system data
load('identified_modal_250.mat');%State-space model
n_x = size(A,1);
n_y = size(C,1);
%Disturbance model dx_w/dt=S*x_w+E*h(t), w=Cw*x_w
load('S.mat');
S = blkdiag(S,zeros(1,1));
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
%% Augmented error dynamics of ESO
n_t = n_x+n_xw;
A0 = zeros(n_t,n_t);
E0 = zeros(n_t,n_h);
C0 = zeros(n_y,n_t);
%
A0(1:n_x,:) = [A H*Cw];
A0(n_x+1:n_x+n_xw,n_x+1:n_x+n_xw) = S;
E0(n_x+1:n_t,:) = E;
C0(:,1:n_x) = C;
Cz = [zeros(n_w,n_x) Cw];%Performance output
%% Sector boundeds of the gain function: phi \in [a,b]
a = 0.7;
b = 1;
%% LMI optimization
c1 = 1e+04;
c2 = 1e+04;
threshold = 1e-06;%Threshold for stric inequalities
cvx_begin sdp
variable  P(n_t,n_t) symmetric
variables Y(n_t,n_y) gam(1,1) tau(1,1)
%
psi11 = P*A0+A0'*P+(Cz'*Cz)-tau*a*b*(C0'*C0);
psi12 = -Y+(tau/2)*(a+b)*C0';
psi13 = P*E0;
psi22 = -tau*eye(n_y);
psi23 = zeros(n_y,n_h);
psi33 = -gam*eye(n_h);
Psi = [psi11  psi12  psi13;
     psi12' psi22  psi23;
     psi13' psi23' psi33];
%
minimize(gam)
subject to
gam>=threshold;
tau>=threshold;
P-c1^(-1)*eye(n_t)>=0;
[-c2^2*eye(n_y) Y';
  Y  -eye(n_t)]<=0;
Psi<=-threshold*eye(size(Psi));
cvx_end
%% Results
gam_opt = gam^0.5;
L = P\Y;
save('ESO.mat','L','gam_opt');