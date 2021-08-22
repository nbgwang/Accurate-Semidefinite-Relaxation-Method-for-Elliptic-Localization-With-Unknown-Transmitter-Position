function [Y,y]=SDP(r,d,S_tilde1,Q1,Q_s1)
% This function realizes the SDP algorithm for jointly
% estimating the unknown object and transmitter positions in the presence
% of receiver position errors using multiple transmitter
%
% Input parameter list:
% S_tilde1:  (K x N), receiver position matrix, N is the number of receivers.         
% r:     (M*N x 1), noisy indirect range measurements.
% d:     (M*N x 1), noisy direct range measurements.
% Q1     (2*M*N x 2*M*N), Covariance matrix after narrow scene
% Q_s:   (K*N x K*N), Covairance matrix of the receiver position errors.
% I:     iterations
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
M=length(d)/N;
      
for i=1:M
    for j=1:N
       b_r(N*(i-1)+j)=0.5*(r(N*(i-1)+j)^2-norm(S_tilde1(:,j))^2);
       b_d(N*(i-1)+j)=0.5*(d(N*(i-1)+j)^2-norm(S_tilde1(:,j))^2);
    end
end
b=[b_r,b_d]';
A_r=[kron(-ones(M,1),S_tilde1'),zeros(M*N,M*K),1/2*kron(ones(M,1),ones(N,1)),...
    zeros(M*N,M),blkdiag(r(1:N),r(N+1:2*N),r(2*N+1:3*N)),kron(-1/2*eye(M),ones(N,1))];
A_d=[zeros(M*N,K),kron(-eye(M),S_tilde1'),zeros(M*N,1),-kron(-1/2*eye(M),ones(N,1)),zeros(M*N,2*M)];
A=[A_r',A_d']';
Sigma=Q1;
Omega=[A'*inv(Sigma)*A,-A'*inv(Sigma)*b;-b'*inv(Sigma)*A,b'*inv(Sigma)*b];
for i=1:2
cvx_clear
cvx_begin sdp
cvx_solver sedumi
cvx_precision best
cvx_quiet(1)
variable y(4*K+10,1)
variable Y(4*K+10,4*K+10) symmetric  
minimize (trace([Y,y;y',1]*Omega))
subject to
[Y,y;y',1]>=0;
trace(Y(1:K,1:K))==y(4*K+1);
for i=1:M
    y(4*K+1+i)==trace(Y(i*K+1:(i+1)*K,i*K+1:(i+1)*K));
    Y(4*K+4+i,4*K+4+i)==trace(Y(1:K,1:K))-2*trace(Y(1:K,i*K+1:(i+1)*K))+trace(Y(i*K+1:(i+1)*K,i*K+1:(i+1)*K));
    y(4*K+7+i)==Y(4*K+4+i,4*K+4+i);
    norm(y(1:K)-y(i*K+1:(i+1)*K))<=y((M+1)*K+M+1+i);
    norm(Y(1:K,(M+1)*K+M+1+i)-Y(i*K+1:(i+1)*K,(M+1)*K+M+1+i))<=Y((M+1)*K+M+1+i,(M+1)*K+M+1+i);
end
cvx_end
u_hat=y(1:K);
T_hat=reshape(y(K+1:(M+1)*K),K,M);
for j=1:N
    Rs(j)=norm(u_hat-S_tilde1(:,j));
    D_h(j,(j-1)*K+1:j*K)=(u_hat-S_tilde1(:,j))';
end
B_h=diag(Rs);  
B_r=kron(eye(M),B_h);
D_r=kron(ones(M,1),D_h);
for i=1:M
    for j=1:N
        Bd(N*(i-1)+j)=norm(T_hat(:,i)-S_tilde1(:,j));
        D_d(N*(i-1)+j,(j-1)*K+1:j*K)=(T_hat(:,i)-S_tilde1(:,j));
    end
end
B_d=diag(Bd);
B=blkdiag(B_r,B_d);
D=[D_r',D_d']';
Sigma=B*Q1*B'+D*Q_s1*D';
Omega=[A'*inv(Sigma)*A,-A'*inv(Sigma)*b;-b'*inv(Sigma)*A,b'*inv(Sigma)*b];
end
