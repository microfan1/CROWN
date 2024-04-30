function [R,mu,Sigma,inv_Sigma,e,f]=Data_Generation(p,N,flag,b,bf,mu_f,cov_f)

k=3; %default factor number

%generate factor return
f=zeros(3,N); %first row: f_t=0
%cov_u = [1-(bf(1,1)^2),0,0;0,1-(bf(2,2)^2),0;0,0,1-(bf(3,3)^2)];
cov_u = cov_f -bf*cov_f*bf';

u = mvnrnd(zeros(3,1),cov_u,N)';
for t=2:N
    f(:,t)=mu_f+bf*f(:,t-1)+u(:,t);
end

%Toeplitz covariance of residuals
s=0.25;
temp=repmat([1:p],p,1);
temp=abs(temp-temp');
Sigma_e=s.^temp; %pxp
e=mvnrnd(zeros(1,p),Sigma_e,N)'; %pxN

temp=zeros(p,1);temp(1)=1+s^2;temp(2)=-s;
temp=toeplitz(temp,temp);temp(1,1)=1;temp(p,p)=1;
inv_V=1/(1-s^2)*temp;


R=(b*f+e)';
Sigma=b*cov_f*b'+Sigma_e;
mu=b*mu_f;
e=e';


%inv_Sigma by Sherman-Morrison-Woodbury Formula
if flag==3
    inv_Sigma=inv_V-inv_V*b/(eye(k)+b'*(inv_V+inv_V')/2*b)*b'*inv_V;
else
    inv_Sigma=inv(Sigma);
end
end