function [X_opt] = joint_waveform_and_filter_design(X_dir,A_theta_t,A_theta_k,P_max)
lambda = 0.1;
alpha = 0.05;

x_0 = X_dir(:);
X_0 = reshape(x_0,T,N);
x = x_0;
for iter=1:100
    %% Update w
    A = 0;
    for i = 1:length(target_DoA)
        A = A + sigma_t(i)*A_theta_t(:,:,i)'*x*x'*A_theta_t(:,:,i);
    end
    B = sigma_u*eye(size(A,1));
    for i = 1:length(interference_DoA)
        B = B + sigma_k(i)*A_theta_k(:,:,i)'*x*x'*A_theta_k(:,:,i);
    end
    w = (pinv(B)*A*x)/(x'*A'*pinv(B)*A*x);

    %% Update x
    U = 0;
    V = 0;
    for kk=1:length(interference_DoA)
        U = U + sigma_k(kk)*(A_theta_k(:,:,kk)'*w*w'*A_theta_k(:,:,kk));
    end

    for tt=1:length(target_DoA)
        V = V + sigma_t(tt)*(A_theta_t(:,:,tt)'*w*w'*A_theta_t(:,:,tt));
    end

    daoham_sinr = 2*(U*x*(x')*V*x - (x'*U*x+sigma_u.*w'*w)*V*x)/(x'*V*x)^2;

    daoham_x = (1-rho)*daoham_sinr + 2*rho*H_tilde'*(H_tilde*x - s) + 2*(1-rho)*lambda*(x-x_0);

    d = x - alpha*daoham_x;

    x(ii) = sqrt(P_max/(N*T))*d(ii)/abs(d(ii));
end

X_opt = reshape(x,[T N]);

end