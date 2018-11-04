%%
%This is a module for TM wave
clear
clc
syms n1 n2 n3 kapa gama delta k d theta lambda omiga x
lam0=300e-6;
n1_0=3.5;
n2_0=1;
n3_0=1;
d0=200e-6;
grating_period=300e-6;
k=2*pi/lambda;
%data_path='data/d0=250_period=350_S11_duty=1:10.dat';
%data_path='data/S11_radial.dat'
data_path='data/S11_duty=1:5.dat';
kapa=n1*k*sin(theta);
gama=sqrt((n1^2-n2^2)*k^2-kapa^2);
delta=sqrt((n1^2-n3^2)*k^2-kapa^2);
f(theta)=n1^2*kapa*d*(n3^2*gama+n2^2*delta)*d/(n2^2*n3^2*kapa^2*d^2-n1^4*gama*delta*d^2);
f(theta)=subs(f,[n1 n2 n3 d],[n1_0 n2_0 n3_0 d0]);
f1(theta)=subs(f,lambda,lam0);
% ezplot(f1(theta),[0 pi/2])
% hold on
% ezplot(tan(sin(theta)*d0*n1_0*2*pi/lam0),[0,pi/2])
% hold off
f2(x)=tan(sin(theta)*d0*n1_0*2*pi/lambda)-f;
f3=sprintf('%s',f2);
cal=strrep(f3,'*','.*');
cal=strrep(cal,'/','./');
cal=strrep(cal,'^','.^');
%%
%pl=parpool(16);
data=[];

% theta=0.00001:0.00001:pi/2;
theta=linspace(0,pi/2,10000);
lambda=900e-6;
n1_0=3.5;
%len=(3000-300)/0.2;
len=(1-0.0001)/0.0001;
a=[];
tic
%for lambda=300e-6:1e-6:900e-6
for i=1:len
    %lambda=300e-6+i*0.2e-6;
    freq=(i/10000+0.0001)*1e12;
    lambda=3e8/freq;
    %fprintf(1, repmat('\b',1,count));
    %delete line before  
    %count=fprintf(1,'current line is : %f',((lambda-300e-6).*100./600e-6));
    %delta=- tan((7.*pi.*sin(theta))./(5000.*lambda)) - (343.*pi.*sin(theta).*((45.*pi.^2)./(4.*lambda.^2) - (49.*pi.^2.*sin(theta).^2)./(4.*lambda.^2)).^(1./2))./(25000000.*lambda.*((21609.*pi.^2)./(80000000.*lambda.^2) - (118433.*pi.^2.*sin(theta).^2)./(400000000.*lambda.^2)));
    delta=eval(cal);
    delta=delta((imag(delta)==0));
    a=[a mean(abs(delta))];
    %beta=(343.*pi.*sin(theta).*(2.*pi.*((49.*sin(theta).^2)./(4.*lambda.^2)).^(1./2) + 2.*((45.*pi.^2)./(4.*lambda.^2) - (49.*pi.^2.*sin(theta).^2)./(4.*lambda.^2)).^(1./2)))./(100000000.*lambda.*((49.*pi.^2.*sin(theta).^2)./(25000000.*lambda.^2) - (2401.*pi.*((45.*pi.^2)./(4.*lambda.^2) - (49.*pi.^2.*sin(theta).^2)./(4.*lambda.^2)).^(1./2).*((49.*sin(theta).^2)./(4.*lambda.^2)).^(1./2))./100000000))-tan((7.*pi.*sin(theta))./(5000.*lambda));
    %a=eval(subs(f3,1));
    targ_theta=theta(abs(delta)<0.005);  
    targ_diff=diff(targ_theta);
    targ_new=[targ_theta(targ_diff>0.001) targ_theta(end)];
    targ_new=targ_new(2:end);
%     if(~isempty(targ_new))
%     if(abs(targ_new(end)-pi/2)<0.01)
%         targ_new=targ_new(1:end-1);
%     end
%     end
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
    for b=beta_0
        data=[data struct('omiga',omiga_0,'beta',b)];
    end
end

    
toc
%delete(pl);
draw_data2