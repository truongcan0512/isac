clear; close all; clc;
N = 16; % number of antenna
K = 4; % number of user
L = 20; % length of the communication frame
power = 1; % Transmit power (P_T)
N_montecarlo = 1; % Number of simulation
SNRdB = -2:2:12;
%%-------------Radar Parameters-------------------
delta=pi/180;
theta=-pi/2:delta:pi/2;
target_DoA=[-20*pi/180, 20*pi/180];
interference_DoA = [-50*pi/180, -10*pi/180, 40*pi/180];
sigma_k = [0.09, 0.25 0.15];
sigma_u = [0.25, 0.25];
sigma_n = 0;



for tt=1:N
    for jj=1:length(theta)
        a(tt,jj)=exp(1j*pi*(tt-N/2)*sin(theta(jj))); % Sterring vector
    end
end

% steer = @(theta) exp(1j*pi*(0:N-1)'*sin(theta));
% theta_t = target_DoA;
% theta_i = interference_DoA;
% % Steering Vectors
% %a_t = steer(theta_t);
% A_t = zeros(N, length(target_DoA));
% for t = 1:length(target_DoA)
%     A_t(:, t) = steer(theta_t(t));
% end
% A_i = zeros(N, length(interference_DoA));
% for k = 1:length(interference_DoA)
%     A_i(:, k) = steer(theta_i(k));
% end


% Steering vector for the target
A0_tmp = zeros(N, length(target_DoA));
for tt=1:N
    for jj=1:length(target_DoA)
        A0_tmp(tt,jj)=exp(1j*pi*(tt-N/2))*sin(target_DoA(jj));
    end
end
for i = 1:length(target_DoA)
    A0(:,:,i) = kron(eye(L), A0_tmp(:,i)*A0_tmp(:,i)');
end
% Steering vector for the interference sources
Ak_tmp = zeros(N, length(interference_DoA));
for tt=1:N
    for jj=1:length(interference_DoA)
        Ak_tmp(tt,jj)=exp(1j*pi*(tt-N/2))*sin(interference_DoA(jj));
    end
end
for i = 1:length(interference_DoA)
    Ak(:,:,i) = kron(eye(L), Ak_tmp(:,i)*Ak_tmp(:,i)');
end

beam_width=9;
l=ceil((target_DoA+pi/2*ones(1,length(target_DoA)))/(delta)+ones(1,length(target_DoA)));
Pd_theta=zeros(length(theta),1);
for ii=1:length(target_DoA)
    Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1)=ones(beam_width,1);
end

l=ceil((interference_DoA+pi/2*ones(1,length(interference_DoA)))/(delta)+ones(1,length(interference_DoA)));
for ii=1:length(interference_DoA)
    Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1)=-0.1*ones(beam_width,1);
end
R = waveform_design_multibm_covmat_new( Pd_theta,N,a,theta,power); % Desired Hermitian positive semidefinite covariance matrix
F = chol(R)'; % Cholesky decomposition


% Radar waveform design with CVX
% eta = 1.8;                            % PAPR upper bound (e.g., max 2x average)
% cvx_begin sdp
% variable R(N, N) hermitian semidefinite
% maximize(sum(diag(A_t'*R*A_t)) - sum(diag(A_i'*R*A_i)))  % total gain at targets
% subject to
% trace(R) <= power;                             % total power constraint
% for n = 1:N
%     R(n,n) <= eta * (1/N) * trace(R);      % PAPR constraint per antenna
% end
% cvx_end
% R = waveform_design_multibm_covmat(Pd_theta,N,L,a,theta,power); % Solve Eq. 10
% F = chol(R)'; % Eq. 11
rho = 1; % Trade-off
amp = sqrt(power);



for nn = 1:N_montecarlo
    H = (randn(K,N)+1j*randn(K,N))/sqrt(2);
    H_tilde = kron(eye(L), H);
    N_pbits = 2*K*L;
    msg_bits = randi([0,1],1,N_pbits);
    Y = reshape(QPSK_mapper(msg_bits),[K,L]); % Received signal
    X_dir = Arbitrary_Com_Rad(H,Y,power,F);   % Eq. 15, Directional Beampattern Design (Dir Strict)
    X_trdoff2 = tradeoff_comrad(rho,H,Y,power,X_dir);           %Trade-off ComRad Total Power
    X_trdoff4 = tradeoff_comrad_per_ant(rho,H,Y,power,X_dir);   %Trade-off ComRad Per Ant
    x_init = X_dir(:);
    y = Y(:);

    %% Main loop
    iter = 100;
    cost = zeros(1, iter);
    alpha = 0.05; % descent coefficient
    x = sqrt(power/(L*N))*ones(L*N,1);
    lambda = 1;
    for it = 1:iter
        %% Update w
        % A = 0;
        % for i = 1:length(target_DoA)
        %     %A = A + sigma_t(i)*A_theta_t(:,:,i)'*x*x'*A_theta_t(:,:,i);
        %     A = A + sigma_u(i)*A0(:,:,i);
        % end
        % B = sigma_n.*eye(size(A,1));
        % 
        % for i = 1:length(interference_DoA)
        %     %B = B + sigma_k(i)*A_theta_k(:,:,i)'*x*x'*A_theta_k(:,:,i);
        %     B = B + sigma_k(i)*Ak(:,:,i)*x*x'*Ak(:,:,i)';
        % end
        % 
        % w = (pinv(B)*A*x)/(x'*A'*pinv(B)*A*x);
        % %% Update w
        A = 0;
        for i = 1:length(target_DoA)
            A = A + sigma_u(i)*A0(:,:,i)'*x*x'*A0(:,:,i);
        end
        B = sigma_n.*eye(size(A,1));
        for i = 1:length(interference_DoA)
            B = B + sigma_k(i)*Ak(:,:,i)'*x*x'*Ak(:,:,i);
        end
        
        w = (pinv(B)*A*x)/(x'*A'*pinv(B)*A*x);

        %% Update x
        % U = 0; V = 0;
        % for i = 1:length(target_DoA)
        %     U = U + sigma_u(i)*A0(:,:,i);
        % end
        % for i = 1:length(interference_DoA)
        %     V = V + sigma_k(i)*Ak(:,:,i);
        % end
        % term = 1/(w'*U*x*x'*U'*w)^2*(w'*U*x*x'*U'*w*V'*w*w'*x - (w'*V*x*x'*V'*w + sigma_n*w'*w)*U'*w*w'*U*x);
        % daoham = 2*(1-rho)*term + 2*rho*H_tilde'*(H_tilde*x - y) ;%+ 2*(1-rho)*lambda*(x - x_init);
        %% Update x
        Q = 0;
        P = 0;
        for kk=1:length(interference_DoA)
            Q = Q + sigma_k(kk)*(Ak(:,:,kk)'*w*w'*Ak(:,:,kk));
        end

        for tt=1:length(target_DoA)
            P = P + sigma_u(tt)*(A0(:,:,tt)'*w*w'*A0(:,:,tt));
        end

        daoham_sinr = 2*(Q*x*(x'*P*x) - (x'*Q*x+sigma_n*w'*w)*(P*x))/((x'*P*x)^2);

        daoham = (1-rho)*daoham_sinr ;%+ 2*rho*H_tilde'*(H_tilde*x - y) + 2*(1-rho)*lambda*(x-x_init);


        x = x - alpha*daoham;
        % for n = 1:N*L
        %     x(n) = sqrt(power/(N*L))*x(n)/abs(x(n));
        %     % if x(n) == 0
        %     %     x(n) = sqrt(power/(N*L));
        %     % else
        %     %     x(n) = sqrt(power/(N*L))*x(n)/abs(x(n));
        %     % end
        % end

        % tu = sigma_u(1)*abs(w'*A*x)^2;
        % tmp = 0;
        % for i = 1:length(interference_DoA)
        %     tmp = tmp + sigma_k(i)*Ak(:,:,i)'*x*x'*Ak(:,:,i);
        % end
        % mau = w'*tmp*w + sigma_n*w'*w;
        mau = 0;
        for i=1:length(target_DoA)
            mau = mau + sigma_u(i)*abs(w'*A0(:,:,i)*x)^2;
        end

        tu = 0;
        for i=1:length(interference_DoA)
            tu = tu + sigma_k(i)*abs(w'*Ak(:,:,i)*x)^2;
        end
        % rho*norm(H_tilde*x - y, 2)^2 + (1-rho)*lambda*norm(x - x_init, 2)^2
        cost(it) =  (1-rho)*tu/mau;
    end
    X_opt = reshape(x, [N L]);
    for ii = 1:length(SNRdB)
        N0 = power/(10^(SNRdB(ii)/10));
        for mm = 1:L
            MUI_arbi(:,mm) = abs(H*X_dir(:,mm)-amp*Y(:,mm)).^2;
            MUI_our(:,mm) = abs(H*X_opt(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff2(:,mm) = abs(H*X_trdoff2(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff4(:,mm) = abs(H*X_trdoff4(:,mm)-amp*Y(:,mm)).^2;
        end
        EMUI_arbi = mean(MUI_arbi,2);
        EMUI_our = mean(MUI_our,2);
        EMUI_trdoff2 = mean(MUI_trdoff2,2);
        EMUI_trdoff4 = mean(MUI_trdoff4,2);

        sumrate_arbi(ii,nn) = sum(log2(1+power./(EMUI_arbi+N0*ones(K,1))));
        sumrate_our(ii,nn) = sum(log2(1+power./(EMUI_our+N0*ones(K,1))));
        sumrate_trdoff2(ii,nn) = sum(log2(1+power./(EMUI_trdoff2+N0*ones(K,1))));
        sumrate_trdoff4(ii,nn) = sum(log2(1+power./(EMUI_trdoff4+N0*ones(K,1))));
        sumrate_lim(ii,nn) = sum(log2(1+power./(N0*ones(K,1))));
    end
end
figure(1)
plot(SNRdB,mean(sumrate_arbi,2),'o-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(SNRdB,mean(sumrate_our,2),'^-','LineWidth',1.5,'MarkerSize',8);
plot(SNRdB,mean(sumrate_trdoff2,2),'*-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_trdoff4,2),'+--','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_lim,2),'v--','LineWidth',1.5,'MarkerSize',8);
legend('Arbi','[9] with multi-target', 'Liu2018 Total', 'Liu2018 Per Ant', 'Zero MUI')
xlabel('Transmit SNR')
ylabel('Communication rate (bps/Hz)')
% legend('Reference Waveform','[9] with multi-target', 'Zero MUI')

figure(2)
plot(theta*180/pi,10*log10(diag(a'*X_opt*X_opt'*a)/real(trace(X_opt*X_opt'))),'LineWidth',1.5);
hold on; grid on;
plot(theta*180/pi,10*log10(diag(a'*X_dir*X_dir'*a)/real(trace(X_dir*X_dir'))),'LineWidth',1.5);
%plot(theta*180/pi,10*log10(diag(a'*X_trdoff2*X_trdoff2'*a)/real(trace(X_trdoff2*X_trdoff2'))),'LineWidth',1.5);
%plot(theta*180/pi,10*log10(diag(a'*X_trdoff4*X_trdoff4'*a)/real(trace(X_trdoff4*X_trdoff4'))),'LineWidth',1.5);
xline(target_DoA*180/pi, 'b-.', 'Linewidth', 1);
xline(interference_DoA*180/pi, 'k-.', 'Linewidth', 1);
legend('[9] with multi-target', 'Reference', 'Liu 2018, Dir, Total Power', 'Liu 2018, Dir, Per Ant')
% legend('[9] with multi-target', 'Reference Waveform')
xlabel('Angle (Degree)')
ylabel('Gain (dB)')


figure(3)
plot(abs(cost))
xlabel('iter')
ylabel('Objective Function Value')