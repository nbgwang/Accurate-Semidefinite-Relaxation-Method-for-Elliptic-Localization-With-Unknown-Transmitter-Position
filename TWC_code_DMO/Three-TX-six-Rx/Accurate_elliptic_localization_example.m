clear all; clc;warning('off'); close all;tic % program initialization.
X=1000;    % This paper has numerical problems, so the scene is reduced by a factor of a thousand
L=1000;    % Number of ensemble runs.           
S=[-2108.1  -2950.6  -401.73   909.87  -3062.2   -116.99;% S and To are randomly generated  
    -3219.8  2036.5  -2405.3   -1233.8  -2553.7  2402.2; % receiver and transmitter positions, respectively
    1819.5   1696.2   2473     2445.3   1387.7   1283.8];
To=[-1115.6  1352.5  2430.5;
    -764.85  3657.7  2385.5;
    1451.6   1863.3  2990.5];
uo=[-1000,500,1500]';                 % Object position
[K,N] = size(S);                      % N is number of receivers and K=dimension
M=3;                                  % Number of transmitters
deltas=1;                             % Receiver position noise power                                
J=diag([5,5,5,40,40,40,20,20,20,10,10,10,30,30,30,15,15,15]);            
Q_s=J*deltas;                         % Covariance matrix of receiver position noise
Q_s1=(1/X)^2*Q_s;                     % Covariance matrix of receiver position noise after narrow scene
a=-10:5:40;                           % Measurements noise power in log-scale

randn('seed',123);                         
noise_rd=randn(2*M*N,L);              % Generating measurement noise

randn('seed',456);
noise_s=randn(K*N,L);                 % Generating receiver position error noise

for iii=1:length(a)                   % Different noise
    fprintf('nsePwr: %d dB\n',a(iii));
    sigma=10^(a(iii)/10);             % Measurement noise power 
   
    for i=1:M
          for j=1:1:N
             r_o(N*(i-1)+j)=norm(uo-To(:,i))+norm(uo-S(:,j)); % true indirect ranges
             d_o(N*(i-1)+j)=norm(To(:,i)-S(:,j));             % true direct ranges
          end
    end
%==========================Construct the covariance matrix=======================
    for i=1:M
        for j=1:1:N
             q_r(N*(i-1)+j)=(r_o(N*(i-1)+j)/sqrt(10^(-0.3)))^2; 
             q_d(N*(i-1)+j)=(d_o(N*(i-1)+j))^2;
        end
    end
    Q_r_bar=diag(q_r);
    Q_d_bat=diag(q_d);
    Qo=blkdiag(Q_r_bar,Q_d_bat);
    qo=sqrt(diag(Qo));
    Q_tilde=0.5*(Qo+qo*qo');
    k_q=2*M*N*sigma/trace(Q_tilde);
    Q=k_q*Q_tilde;                   % Covariance matrix of indirect and direct range noise
    Q1=(1/X)^2*Q;                    % Covariance matrix after narrow scene
      
    sum_rank=0;sum_u=0;              % Initialize
    sum_crb_u=0; 
    
%==========================Generates zero mean noise==========================
    noise_rd1=noise_rd-mean(noise_rd,2)*ones(1,L);
    error_rd=sqrtm(Q)*noise_rd1;
    error_r=reshape(error_rd(1:M*N,:),M*N,1,L);
    error_d=reshape(error_rd(M*N+1:2*M*N,:),M*N,1,L);

    noise_s1=noise_s-mean(noise_s,2)*ones(1,L);
    error_s=sqrtm(Q_s)*noise_s1;
    error_s=reshape(error_s,K,N,L);
    
parfor ii=1:L                                           % Monte Carlo number
     S_tilde=S+error_s(:,:,ii);                         % Noisy receiver positions
     S_tilde1=1/X*S_tilde;                              % Receiver position with noise narrow scene
   
     r=1/X*(r_o'+error_r(:,:,ii));                       % Indirect range measurement after narrow scene
     d=1/X*(d_o'+error_d(:,:,ii));                       % Direct range measurement after narrow scene
     
     [Y,y]=SDP(r,d,S_tilde1,Q1,Q_s1);                   % SDP solution  
     landa=eig(Y(1:4*K,1:4*K));                         % Eigenvalue decomposition for Y(1:4*K,1:4*K)
     landa=sort(landa,'descend');                       % The eigenvalues are sorted from large to small
     if (landa(1)/landa(2))>10^5                                   
          rank_sdp=1;
          sum_rank=sum_rank+rank_sdp;                   % Count the number of rank1
          rank_ut(iii,ii)=1;
     else
          rank_sdp=0;
          rank_ut(iii,ii)=0;
     end
     u_sdp_u(:,ii)=X*y(1:K);                            % Object position estimation
     sum_u=sum_u+(norm(u_sdp_u(:,ii)-uo))^2;          
end
sumu_rank_total(iii)=sum_rank;              
error_total_sdp_u(iii)=sum_u;                         

MSE=CRLB(To,uo,S,Q,Q_s);                                % CRLB evaluation
y_crb_u(iii)=10*log10(trace(MSE(1:K,1:K)));
end

MSE_average_u=1/L*error_total_sdp_u;                    % average MSE of object position
y_sdp_u=10*log10(MSE_average_u);                                   


figure(1);
plot(a,y_sdp_u,'r^');
hold  on;
plot(a,y_crb_u,'-k','MarkerSize',8);
hold off;
grid on; xlabel('10lg(\sigma^2(m^2))'); ylabel('10lg(MSE(m^2))');
legend('SDP','CRLB')
toc









