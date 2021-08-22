function MSE=CRLB(To,uo,S,Q,Q_s)
% This function computes the CRLB of the joint estimation of the unknown object
% and transmitter positions in the prsence of receiver position errors using multiple transmitters.
%
% Input parameter list:
% S:    (K x N), receiver position matrix, N is the number of receivers.         
% To:   (K x M), transmitter position, M is the number of transmitters. 
% uo:   (K x 1), object position.
% Q     (2*M*N x 2*M*N), Covariance matrix after narrow scene
% Q_s:  (K*N x K*N), Covairance matrix of the receiver position errors.
%
% Output parameter list:
% CRLB:  (K*(M+N+1) x K*(M+N+1)), CRLB matrix of the object and transmitter
%        and receivers positions estimate
%            
% Reference:
% Ruichao Zheng, Gang Wang and K. C. Ho, "Accurate Semidefinite Relaxation Method
% for Elliptic Localization With Unknown Transmitter Position," 
% IEEE Transactions on Wireless Communications
% vol. 20, no. 4, pp. 2746-2760, April 2021.
% R. Zheng, G. Wang and K. C. Ho

[K,N]=size(S);
M=3;
for i=1:M
      for j=1:N
          deltaru((i-1)*N+j,:)=(uo-S(:,j))'/norm(uo-S(:,j))-(To(:,i)-uo)'/norm(To(:,i)-uo);
          delta((i-1)*N+j,:)=(To(:,i)-S(:,j))'/norm(To(:,i)-S(:,j));
          deltart(N*(i-1)+j,(i-1)*K+1:i*K)=(To(:,i)-uo)'/norm(To(:,i)-uo);
          deltars(N*(i-1)+j,(j-1)*K+1:j*K)=(S(:,j)-uo)'/norm(S(:,j)-uo);
          deltads(N*(i-1)+j,(j-1)*K+1:j*K)=(S(:,j)-To(:,i))'/norm(S(:,j)-To(:,i));
      end
end
deltadt=blkdiag(delta(1:N,1:K),delta(N+1:2*N,1:K),delta(2*N+1:3*N,1:K));
A=[zeros((M+1)*K,(M+N+1)*K);zeros(N*K,(M+1)*K),inv(Q_s)];
deltard=[deltaru,deltart,deltars;zeros(M*N,K),deltadt,deltads];
MSE=inv(deltard'*inv(Q)*deltard+A);
end