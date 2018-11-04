%%
%This is a module for TM wave
% clear
% clc
% syms n1 n2 n3 kapa gama delta k d theta lambda omiga 
% lam0=900e-6;
% n1_0=3.5;
% n2_0=1;
% n3_0=1;
% d0=200e-6;
% k=2*pi/lambda;
% kapa=n1*k*sin(theta);
% gama=sqrt((n1^2-n2^2)*k^2-kapa^2);
% delta=sqrt(kapa^2);
% f(theta)=n1^2*kapa.*d*(n3^2*gama+n2^2*delta)*d/(n2^2*n3^2*kapa^2*d^2-n1^4*gama*delta*d^2);
% f(theta)=subs(f,[n1 n2 n3 d lambda],[n1_0 n2_0 n3_0 d0 lam0]);
% ezplot(f(theta),[0 pi/2])
% hold on
% ezplot(tan(sin(theta).*200e-6.*3.5.*2.*pi./lam0),[0,pi/2])
% hold off
% f2=f-tan(sin(theta).*200e-6.*3.5.*2.*pi./lambda);
%%
clear
clc
pl=parpool(16);
data=[];
count=0;
theta=0.00001:0.00001:pi/2;
lambda=900e-6;
n1_0=3.5;
len=(900-300)/0.2;
a=[];
tic
%for lambda=300e-6:10e-6:900e-6
parfor i=1:len
    lambda=300e-6+i*0.2e-6;
    %fprintf(1, repmat('\b',1,count));
    %delete line before  
    %count=fprintf(1,'current line is : %f',((lambda-300e-6).*100./600e-6));
    delta=- tan((7.*pi.*sin(theta))./(5000.*lambda)) - (7.*pi.*sin(theta).*((45.*pi.^2)./(4.*lambda.^2) - (49.*pi.^2.*sin(theta).^2)./(4.*lambda.^2)).^(1./2))./(6250000.*lambda.*((9.*pi.^2)./(5000000.*lambda.^2) - (49.*pi.^2.*sin(theta).^2)./(12500000.*lambda.^2)));
    a=[a mean(delta)];
    targ_theta=theta(abs(delta)<0.01);  
    targ_diff=diff(targ_theta);
    targ_new=[targ_theta(targ_diff>0.001) targ_theta(end)];
    %delete useless points
%     len_targ=length(targ_theta);
%     index_targ=1;
%     targ_new=[];
%     while(index_targ<=len_targ)
%         targ_new=[targ_new mean(targ_theta(abs(targ_theta-targ_theta(index_targ))<0.01))];
%         index_targ=index_targ+length(targ_theta(abs(targ_theta-targ_theta(index_targ))<0.01));
%     end
    omiga_0=2.*pi.*3e8./lambda;
    beta_0=n1_0.*2.*pi./lambda.*cos(targ_new);
    data=[data struct('omiga',omiga_0,'beta',beta_0)];
end
toc
delete(pl)