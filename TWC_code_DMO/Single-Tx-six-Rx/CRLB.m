function MSE=CRLB(to,uo,S,Q,Q_s)
% This function computes the CRLB of the joint estimation of the unknown object
% and transmitter positions in the prsence of receiver position errors using single transmitter.
%
% Input parameter list:
% S:    (K x N), receiver position matrix, N is the number of receivers.         
% to:   (K x 1), transmitter position, M is the number of transmitters. 
% uo:   (K x 1), object position.
% Q     (2*M*N x 2*M*N), Covariance matrix after narrow scene
% Q_s:  (K*N x K*N), Covairance matrix of the receiver position errors.
%
% Output parameter list:
% CRLB:  (K*(N+2) x (K*(N+2), CRLB matrix of the object and transmitter
%        and receivers positions estimate
%            
% Reference:
% Ruichao Zheng, Gang Wang and K. C. Ho, "Accurate Semidefinite Relaxation Method
% for Elliptic Localization With Unknown Transmitter Position," 
% IEEE Transactions on Wireless Communications
% vol. 20, no. 4, pp. 2746-2760, April 2021.
% R. Zheng, G. Wang and K. C. Ho

[K,N]=size(S);
p=(to-uo)'/norm(to-uo);
datart=kron(ones(N,1),p);
for j=1:N
    datau(j,:)=(uo-S(:,j))'/norm(uo-S(:,j));
    datadt(j,:)=(to-S(:,j))'/norm(to-S(:,j));
    datars(j,(j-1)*K+1:j*K)=(S(:,j)-uo)'/norm(S(:,j)-uo);
    datads(j,(j-1)*K+1:j*K)=(S(:,j)-to)'/norm(S(:,j)-to);
end
dataru=datau-datart;
A=[zeros(2*K,(N+2)*K);zeros(N*K,2*K),inv(Q_s)];
datard=[dataru,datart,datars;zeros(N,K),datadt,datads];
MSE=inv(datard'*inv(Q)*datard+A);
end