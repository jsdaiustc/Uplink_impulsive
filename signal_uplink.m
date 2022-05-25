function [Y,S,H]=signal_uplink(K,N,L,T,SNR)

Ps=sqrt((   10.^(SNR/10)   )/2);
X=Ps*( randn(L,K)+1i*randn(L,K)  );   
DOA=(randn(1,L)-0.5)*90;
A = exp(-1i*pi*(0:N-1)'*sind(DOA));
H=A*X;
S=( randn(K,T)+1i*randn(K,T)  );



%%  GMM noise
   Mu = [0;0];
   c=0.1;
   for i=1:T
       Sigma = cat(3,[1],[100]);
       P=[1-c,c];
       gm = gmdistribution(Mu,Sigma,P);
       r1=random(gm,N);
       r2=random(gm,N);
       noise(:,i)=sqrt(1/2)*(r1+1i*r2);
   end
 Y=H*S + noise;





