function [Y,y]=SDP(r,d,S_tilde1,Q1,Q_s1)
% This function realizes the SDP algorithm for jointly
% estimating the unknown object and transmitter positions in the presence
% of receiver position errors using signle transmitter
%
% Input parameter list:
% S_tilde1:  (K x N), receiver position matrix, N is the number of receivers.         
% r:     (N x 1), noisy indirect range measurements.
% d:     (N x 1), noisy direct range measurements.
% Q1     (2*N x 2*N), Covariance matrix after narrow scene
% Q_s:   (K*N x K*N), Covairance matrix of the receiver position errors.
% 
% Output parameter list:
% y:  Estimated object and transmitter positions.
%
% The program can be used for 2D(K=2) or 3D(K=3) localization
%
% Reference:
% Ruichao Zheng, Gang Wang and K. C. Ho, "Accurate Semidefinite Relaxation Method
% for Elliptic Localization With Unknown Transmitter Position," 
% IEEE Transactions on Wireless Communications
% vol. 20, no. 4, pp. 2746-2760, April 2021.
% R. Zheng, G. Wang and K. C. Ho

[K,N]=size(S_tilde1);                           
for j=1:N
   b_r(j)=0.5*(r(j)^2-norm(S_tilde1(:,j))^2);   
   b_d(j)=0.5*(d(j)^2-norm(S_tilde1(:,j))^2);
end
b=[b_r,b_d]';
A=[-S_tilde1',zeros(N,K),0.5*ones(N,1),zeros(N,1),blkdiag(r),-0.5*ones(N,1);
    zeros(N,K),-S_tilde1',zeros(N,1),0.5*ones(N,1),zeros(N,2)];
Sigma=Q1;
Omega=[A'*inv(Sigma)*A,-A'*inv(Sigma)*b;-b'*inv(Sigma)*A,b'*inv(Sigma)*b];
for i=1:2 
cvx_clear
cvx_begin sdp
cvx_solver sedumi
cvx_precision best
cvx_quiet(1)
variable y(2*K+4,1)
variable Y(2*K+4,2*K+4) symmetric  
minimize (trace([Y,y;y',1]*Omega))
subject to
[Y,y;y',1]>=0;
Y(2*K+3,2*K+3)==trace(Y(1:K,1:K))+trace(Y(K+1:2*K,K+1:2*K))-2*trace(Y(1:K,K+1:2*K));
trace(Y(1:K,1:K))==y(2*K+1);
trace(Y(K+1:2*K,K+1:2*K))==y(2*K+2);
trace(Y(2*K+3,2*K+3))==y(2*K+4);
norm(y(1:K)-y(K+1:2*K))<=y(2*K+3);
norm(Y(1:K,2*K+3)-Y(K+1:2*K,2*K+3))<=Y(2*K+3,2*K+3);
cvx_end
u_hat=y(1:K);
t_hat=y(K+1:2*K);
for j=1:N
    Rs(j)=norm(u_hat-S_tilde1(:,j));
    Dd(j)=norm(t_hat-S_tilde1(:,j));
end
B_r=diag(Rs);
B_d=diag(Dd);
B=blkdiag(B_r,B_d);
for j=1:N
    D_r(j,(j-1)*K+1:j*K)=(u_hat-S_tilde1(:,j))';
    D_d(j,(j-1)*K+1:j*K)=(t_hat-S_tilde1(:,j))';
end
D=[D_r',D_d']';
Sigma=B*Q1*B'+D*Q_s1*D';
Omega=[A'*inv(Sigma)*A,-A'*inv(Sigma)*b;-b'*inv(Sigma)*A,b'*inv(Sigma)*b];
end


