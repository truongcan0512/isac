function [x,w,cost] = joint_trade_off_CM( ...
    rho, ...    % trade-off params
    H, ...      % Channel
    Y, ...      % Y=HX+W
    N, ...      % symbol times
    T, ...      % No.transmit antenas
    R, ...      % No.receive antenas
    P_max, ...  % transmit power
    target_DoA, ...     % targets
    interference_DoA, ... % interferences
    sigma_t, ...    % complex amplitudes of targets
    sigma_k, ...    % complex amplitudes of interferencess
    sigma_n, ...    % complex amplitudes of noise
    max_iter, ...   % max iteration
    alpha ...      % descent rate
    )

%Steering vectors of targets
for tt=1:T
    for jj=1:length(target_DoA)
        A_tt(tt,jj)=exp(-1j*pi*(tt-1)*sin(target_DoA(jj)));
    end
end

for rr=1:R
    for jj=1:length(target_DoA)
        A_rt(rr,jj)=exp(-1j*pi*(rr-1)*sin(target_DoA(jj)));
    end
end

for tt=1:length(target_DoA)
    A_theta_t(:,:,tt) = kron(eye(N),A_rt(:,tt)*A_tt(:,tt)');
end


%Steering vectors of interferences
for tt=1:T
    for jj=1:length(interference_DoA)
        A_tk(tt,jj)=exp(-1j*pi*(tt-1)*sin(interference_DoA(jj)));
    end
end

for rr=1:R
    for jj=1:length(interference_DoA)
        A_rk(rr,jj)=exp(-1j*pi*(rr-1)*sin(interference_DoA(jj)));
    end
end

for kk=1:length(interference_DoA)
    A_theta_k(:,:,kk) = kron(eye(N),A_rk(:,kk)*A_tk(:,kk)');
end


H_tilde = kron(eye(N),H);
s = Y(:);

x = sqrt(P_max)*ones(T*N,1); % CM domain
cost = zeros(1, max_iter);
for iter=1:max_iter
    %% Update w
    A = 0;
    for i = 1:length(target_DoA)
        A = A + sigma_t(i)*A_theta_t(:,:,i);
    end
    B = sigma_n.*eye(size(A,1));

    for i = 1:length(interference_DoA)
        B = B + sigma_k(i)*A_theta_k(:,:,i)*x*x'*A_theta_k(:,:,i)';
    end

    w = (pinv(B)*A*x)/(x'*A'*pinv(B)*A*x);

    %% Update x
    Q = 0;
    P = 0;
    for kk=1:length(interference_DoA)
        Q = Q + sigma_k(kk)*(A_theta_k(:,:,kk)'*w*w'*A_theta_k(:,:,kk));
    end

    for tt=1:length(target_DoA)
        P = P + sigma_t(tt)*(A_theta_t(:,:,tt)'*w*w'*A_theta_t(:,:,tt));
    end

    daoham_sinr = 2*(Q*x*x'*P*x - (x'*Q*x+sigma_n*w'*w)*(P*x))/((x'*P*x)^2);
    
    daoham_x = 2*rho*H_tilde'*(H_tilde*x - s)  + (1-rho)*daoham_sinr;% +  2*(1-rho)*lambda*(x-x_0);
    
    % Project Gradient Descent method
    x = x - alpha*daoham_x;
    
    for ii=1:length(x)
        % if abs(x(ii))^2 > P_max
        %     x(ii) = sqrt(P_max)*x(ii)/abs(x(ii)); % Total energy constraint
        % end

        % Constant modulus constraint
        if abs(x(ii)) == 0
            x(ii) = sqrt(P_max/(N*T));
        else
            x(ii) = sqrt(P_max/(N*T))*x(ii)/abs(x(ii));
        end
    end

    %% Tinh cost function

    % Power of Interference
    power_i = 0;
    for i=1:length(interference_DoA)
        power_i = power_i + (sigma_k(i)*abs(w'*A_theta_k(:,:,i)*x))^2;
    end

    power_i = power_i + sigma_n*w'*w;
    % Power of Signal
    power_t = 0;
    for i=1:length(target_DoA)
        power_t = power_t + (sigma_t(i)*abs(w'*A_theta_t(:,:,i)*x))^2;
    end
    cost(iter) =  rho*norm(H_tilde*x - s, 2)^2 + (1-rho)*power_i/power_t;% + (1-rho)*lambda*norm(x - x_0, 2)^2;
end