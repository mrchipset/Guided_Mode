%%
%This is a module for TM wave
clear
clc
syms n1 n2 n3 kapa gama delta k d theta lambda omiga
lam0=300e-6;
n1_0=3.5;
n2_0=1;
n3_0=1;
d0=200e-6;
k=2*pi/lambda;
kapa=n1*k*sin(theta);
gama=sqrt((n1^2-n2^2)*k^2-kapa^2);
delta=sqrt((n1^2-n3^2)*k^2-kapa^2);
f(theta)=n1^2*kapa*d*(n3^2*gama+n2^2*delta)*d/(n2^2*n3^2*kapa^2*d^2-n1^4*gama*delta*d^2);
f(theta)=subs(f,[n1 n2 n3 d],[n1_0 n2_0 n3_0 d0]);
f1(theta)=subs(f,lambda,lam0);
ezplot(f1(theta),[0 pi/2])
hold on
ezplot(tan(sin(theta)*d0*n1_0*2*pi/lam0),[0,pi/2])
hold off
f2(k)=tan(sin(theta)*d0*n1_0*2*pi/lambda)-f;
f3=sprintf('%s',f2);
cal=strrep(f3,'*','.*');
cal=strrep(cal,'/','./');
cal=strrep(cal,'^','.^');

% data=[];
% count=0;
% for lam0=300e-6:0.5e-6:900e-6
%     fprintf(1, repmat('\b',1,count));
%     %delete line before  
%     count=fprintf(1,'current line is : %f',((lam0-300e-6)*100/600e-6));
%     f3=subs(f2,lambda,lam0);
%     b=0.0001:0.0001:pi/2;
%     a=eval(subs(f3,b));
%     %a=eval(subs(f3,1));
%     targ_theta=b(abs(a)<0.01);  
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