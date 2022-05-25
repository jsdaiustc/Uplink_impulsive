function H=offgrid_Bayesian_uplinkImpulsivefast(Y,S,N_point,etc)
hatH=Y/S; 
S_inv=pinv(S);
[N,K]=size(hatH);
[T,~]=size(S_inv);

norm_h=norm(hatH,'fro')/sqrt(N*K);
hatH=hatH/norm_h;


%% initialization
search_area1=[-1:2/N_point:1];   
search_area=ones(K,1)*search_area1;
reslu=search_area1(3)-search_area1(2);   
A=exp(-1i*pi*(0:N-1)'*search_area1)/sqrt(N);

for k=1:K
    phi_e(:,:,k)=kron(S_inv(:,k).',eye(N));
    AK(:,:,k)=A;    
end

a=0.0001;b=0.0001;
maxiter=100;
tol=1e-5;
rho=0.01;
converged = false;
iter = 0;
gamma_inv =ones(length(search_area1),K);
alpha=1;
mu_e=zeros(T*N,1);
Z=mu_e;
delta=ones(T*N,1);
sum_phi_e= kron( conj(S_inv)*S_inv.' ,eye(N));

for kk=1:K
    SVD(:,kk) = svd(phi_e(:,:,kk),'econ');
end
eta=max(max(SVD))^2+0.00001;

% Y_abs=abs(Y(:));
% Y(Y_abs>median(Y_abs)*10)=0;


% algrithm begins
while ~converged
    
    
    %% calculate mu and Sigma q(xk)
    for k=1:K
        Phi_gamma = AK(:,:,k) *  diag(gamma_inv(:,k));
        V_temp= 1/alpha*eye(N) + Phi_gamma * AK(:,:,k)';
        Sigma_x (:,:,k)= diag(gamma_inv(:,k)) -Phi_gamma' * (V_temp \Phi_gamma);
        mu_x(:,k) = alpha * Sigma_x(:,:,k) * AK(:,:,k)' * (hatH(:,k) -phi_e(:,:,k)*mu_e);
    end
    
    %% update q(e)
    for k=1:K
        q(:,k)=phi_e(:,:,k)'*(hatH(:,k)- AK(:,:,k)*mu_x(:,k));
    end
    Sigma_e= 1./(alpha*K*eta+rho*delta);
    mu_e=alpha*Sigma_e.*(K*eta*Z+sum(q,2)-sum_phi_e*Z);
    
    
    %% %%%%%%update alpha
    for k=1:K
        resid(:,k)=hatH(:,k)-AK(:,:,k)*mu_x(:,k)-phi_e(:,:,k)*mu_e;
        norm1(:,k)=norm(resid(:,k), 'fro')^2+eta*sum((mu_e-Z).^2)+eta*sum(Sigma_e)+sum(diag(AK(:,:,k)*Sigma_x(:,:,k)*AK(:,:,k)'))+2*real( (mu_e-Z )'* phi_e(:,:,k)'*  (-hatH(:,k)+ AK(:,:,k)*mu_x(:,k)+phi_e(:,:,k)*Z));
    end
    alpha=( K*N + a)/( b + real( sum(norm1,2)) ) ;
    
    
    %% %%%%update gamma
    gamma_last = gamma_inv;
    for k=1:K
        sum_mu(:,k)=sum( mu_x(:,k).*conj(mu_x(:,k) ), 2);
        temp1= sum_mu(:,k) + real(diag(Sigma_x(:,:,k))) ;
        gamma_inv(:,k)=    ( b+ real(temp1)  )/(a+1);
    end
    
    
    
    %% update delta
    sum_emu=sum( mu_e.*conj(mu_e), 2);
    temp=sum_emu + real(Sigma_e);
    delta=   (a+1)./( b+ rho*real(temp));
    
    
    %% update Z
    Z=mu_e;
    
    
    
    %% %%%%%%%%grid refine
    
    for k=1:K
        Pm(:,k)=sum( mu_x(:,k).*conj(mu_x(:,k)), 2);
        [~,sort_ind]=sort(Pm(:,k), 'descend');
        index_amp = sort_ind(1:etc);
        tempPS=AK(:,:,k)*  Sigma_x(:,index_amp,k) ;
        df=zeros(length(index_amp),1);
        for j=1:length(index_amp)
            ii=index_amp(j);
            ai= exp(-1i*pi*(0:N-1)'*search_area(k,ii))/sqrt(N);
            mut=mu_x(ii,k);
            Sigmat=Sigma_x(:,ii,k);
            c1=mut*mut' +  Sigmat(ii);
            c1=abs(c1)*(-alpha);
            Yti(:,k)=resid(:,k)+  AK(:, ii,k)*mu_x(ii,k);
            c2=  (  tempPS(:,j) - AK(:,ii,k)*Sigmat(ii) )  -Yti(:,k)*(mut');
            c2= c2*(-alpha);
            phii=AK(:,ii,k);
            sinta= search_area(k,ii);
            costa=cos(asin(sinta));
            c3=(-1i*pi* costa)*[0:N-1]';
            tt1=  (c3.*ai);
            f1=2*real(tt1'*phii*c1)  +  2*real( tt1'*c2);
            df(j)=f1;
        end
        ddff=sign(df.')*reslu/100;
        search_area(k,index_amp) = search_area(k,index_amp) +  ddff;
        F_active=exp(-1i*pi*(0:N-1)'* search_area(k,index_amp))/sqrt(N);
        AK(:,index_amp,k)=F_active;
        
    end
    
    
    
    % stopping criteria
    erro=norm(gamma_inv - gamma_last)/norm(gamma_last);
    if erro < tol || iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end


for k=1:K
    H(:,k)=AK(:,:,k)*mu_x(:,k)*norm_h;
end
