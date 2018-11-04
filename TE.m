%%
%Draw define omiga and beta TE
%
clear
clc
syms n1 n2 n3 kapa gama delta k d theta lambda
lam0=350e-6;
n1_0=3.5;
n2_0=1;
n3_0=1;
d_0=200e-6;
%omiga=3e8/lambda*2*pi;
k=2*pi/lambda;
kapa=n1*k*sin(theta);
gama=sqrt((n1^2-n2^2)*k^2-kapa^2);
delta=sqrt((n1^2-n2^2)*k^2-kapa^2);
f(theta)=kapa*d*(gama+delta)*d/(kapa^2*d^2-delta*gama*d^2);
%f(kapa)=subs(f,lambda,2*pi*n1*sin(theta)/kapa);
f2=f-tan(sin(theta)*d_0*n1_0*2*pi/lambda);
f2=subs(f2,[n1 n2 n3 d],[n1_0 n2_0 n3_0 d_0]);
f(theta)=subs(f,[n1 n2 n3 d lambda],[n1_0 n2_0 n3_0 d_0 lam0]);
ezplot(f(theta),[0 pi/2])
hold on
%figure
ezplot(tan(sin(theta)*200e-6*3.5*2*pi/lam0),[0,pi/2])
hold off

%%
%This is a test for lambda expression
%parpool(4);
% data=[];
% count=0;
% b=0.0001:0.0001:pi/2;
% a=zeros(1,length(b));
% omiga_0=0;
% beta_0=0;
% for lam0=300e-6:0.5e-6:900e-6  
%     fprintf(1, repmat('\b',1,count));
%     %delete line before  
%     count=fprintf(1,'current line is : %d',floor(lam0/100e-6)-2);
%     f3=subs(f2,lambda,lam0);     
%     a=eval(subs(f3,b));
%     %a=eval(subs(f3,1));
%     targ_theta=b(abs(a)<0.01);
%     %%
%     %delete useless points
%     len_targ=length(targ_theta);
%     index_targ=1;
%     targ_new=[];
%     while(index_targ<=len_targ)
%         targ_new=[targ_new mean(targ_theta(abs(targ_theta-targ_theta(index_targ))<0.01))];
%         index_targ=index_targ+length(targ_theta(abs(targ_theta-targ_theta(index_targ))<0.01));
%     end
%     omiga_0=2*pi*3e8/lam0;
%     beta_0=n1_0*2*pi/lam0*sin(targ_new);
%     data=[data struct('omiga',omiga_0,'beta',beta_0)];
% end