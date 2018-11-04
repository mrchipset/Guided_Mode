%%
%Dump the data
omiga=[];
beta=[];
p=[];
p_2=[];
p_3=[];
min0=1e12;
min1=1e12;
min_omiga=0;
min_omiga1=0;

c=3e8;
n1=1;
n2=3.5;
len=length(data);
%S11=load('data/d0=250_period=350_S11_duty=1:10.dat');
S11=load(data_path);
% pl=parpool(16);
% parfor i=1:len
%     o_data=data(i);
%     o_beta=sort(o_data.beta);
%     dest_beta=[];
%     o_omiga=[];
%     if(length(o_beta)>3)
%         dest_beta=o_beta(1:3);
%         o_omiga=[o_data.omiga o_data.omiga o_data.omiga]
%     else
%         for j=o_beta
%            dest_beta=[dest_beta j];
%            o_omiga=[o_omiga o_data.omiga]
%         end
%     end
d_beta=[];
d_omiga=[];
for i=data  
        omiga=[omiga i.omiga];
        beta=[beta i.beta];
%         p=[p 2*pi/300e-6];
%         p_2=[p_2 4*pi/300e-6];
%         p_3=[p_3 6*pi/300e-6];
%         if(abs(j-2*pi/300e-6)<min0)
%             min0=abs(j-2*pi/300e-6);
%             min_omiga=i.omiga;
%         end
%         
%         if(abs(j-2*pi/150e-6)<min1)
%             min1=abs(j-2*pi/150e-6);
%             min_omiga1=i.omiga;
%         end
            
end
p_beta=[];
p_omiga=[];
p_data=[];
f=fopen('log/grating_cal.txt','w');
fprintf(f,'Grating Beta,\tCalculation Beta,\tCalculation Freq\n');
for i=[2*pi/grating_period 4*pi/grating_period 6*pi/grating_period]
    d_beta=beta((abs(beta-i)/i)<0.01);
    d_omiga=omiga((abs(beta-i)/i)<0.01);
%     d_beta=d_beta(diff(d_omiga)>0.1*1e12);
%     d_omiga=d_omiga(diff(d_omiga)>0.1*1e12);
    p_beta=[p_beta d_beta];
    p_omiga=[p_omiga d_omiga];
    for j=1:length(d_beta)
        fprintf(f,'%.4f,\t%.4f,\t%.4f\n',i,d_beta(j),d_omiga(j));
        p_data=[p_data;i,d_beta(j),d_omiga(j)];
    end
      
%       d_beta=[beta(abs(beta-i)/i<0.03) beta(abs(beta-i)./beta<0.03)];
%       d_omiga=[omiga(abs(beta-i)/i<0.03) omiga(abs(beta-i)./beta<0.03)];
%       dd_data=[];
%       for j=1:length(d_beta)
%         dd_data=[dd_data;i,d_beta,d_omiga];
%       end
     
end
fclose(f);

p_data=drop_data1(p_data);
p_omiga=p_data(:,3);
p_beta=p_data(:,1);
f=fopen('log/data_dropped_grating_cal.txt','w');
fprintf(f,'Grating Beta,\tCalculation Beta,\tCalculation Freq\n');
for i=1:length(p_data)
    fprintf(f,'%.4f,\t%.4f,\t%.4f\n',p_data(i,1),p_data(i,2),p_data(i,3));
end
fclose(f);


%p_beta=sort(p_beta);
%p_omiga=sort(p_omiga);

% delete(pl);
omiga=omiga/2/pi;
p_omiga=p_omiga/2/pi;
c=c/2/pi;

% figure
% plot(omiga,beta,'.k');
% hold on
% plot(omiga,omiga*n1/c);
% plot(omiga,omiga*n2/c);
% hold off
% xlabel('\omega');
% ylabel('\beta');


figure
plot(beta,omiga,'.k');
hold on

p=ones(100,1);
p_2=ones(100,1);
p_3=ones(100,1);
p=p*2*pi/grating_period;
p_2=p_2*4*pi/grating_period;
p_3=p_3*6*pi/grating_period;
plot(p,linspace(0,max(omiga),length(p)));
plot(p_2,linspace(0,max(omiga),length(p_2)));
plot(p_3,linspace(0,max(omiga),length(p_3)));
plot(beta,beta*c/n1)
plot(beta,beta*c/n2)
plot(p_beta,p_omiga,'*r');

hold off
axis([0 7e4 0 1e12])
ylabel('\omega');
xlabel('\beta');


freq=S11(:,2);
ref=S11(:,3);
tar_ref=ref(diff(sign(diff(ref)))~=0);
tar_freq=freq(diff(sign(diff(ref)))~=0);
len=length(tar_freq);
p_omiga=p_omiga/1e12;
tar_freq=tar_freq/1e12;
figure
plot(freq,ref,'--k');
hold on
res=[];
f=fopen('log/cal_sim_point.txt','w');
fprintf(f,'Calculation Freq,\tSimulation Freq,\tSimulation Reflection\n');
c_data=[];
for i=p_omiga'
    for index=1:len
        if(abs(tar_freq(index)-i)/i<0.05||abs(tar_freq(index)-i)/tar_freq(index)<0.05)
%             fprintf(f,'%.4f,\t%.4f,\t%.4f\n',i,tar_freq(index),tar_ref(index));
%             plot(tar_freq(index)*1e12,tar_ref(index),'*r');
            c_data=[c_data;i,tar_freq(index),tar_ref(index)];
        end
    end
end

c_data=drop_data2(c_data);
plot(c_data(:,2)*1e12,c_data(:,3),'*r');
for i = 1:length(c_data)
    fprintf(f,'%.4f,\t%.4f,\t%.4f\n',c_data(i,1),c_data(i,2),c_data(i,3));
end
fclose(f);
hold off


%%Drop complex data point
%
%
%

function res=drop_data1(data)
    o_data=[];
    tmp=data;
    d_data=diff(tmp(:,3));
    min_index=0;
    delta=99999;
    for i=1:length(d_data)
        if(d_data(i)>0.01e12)
            if(min_index~=0)
                if(delta>abs(data(i,1)-data(i,2)))
                    min_index=i;
                end
                o_data=[o_data;data(min_index,:)];
                min_index=0;
                delta=99999;
            else
                o_data=[o_data; data(i,:)];
            end
        else
            if(delta>abs(data(i,1)-data(i,2)))
                delta=abs(data(i,1)-data(i,2));
                min_index=i;
            end
        end
        
    end
    
    if(delta>abs(data(i+1,1)-data(i+1,2)))
        o_data=[o_data;data(i+1,:)];
    else
        o_data=[o_data;data(min_index,:)];
    end
    
    res=o_data;
end





function res=drop_data2(data)
    o_data=[];
    tmp=data;
    d_data=diff(tmp);
    min_index=0;
    delta=99999;
    for i=1:length(d_data)
        if(d_data(i)~=0)
            if(min_index~=0)
                if(delta>abs(data(i,1)-data(i,2)))
                    min_index=i;
                end
                o_data=[o_data;data(min_index,:)];
                min_index=0;
                delta=99999;
            else
                o_data=[o_data; data(i,:)];
            end
        else
            if(delta>abs(data(i,1)-data(i,2)))
                delta=abs(data(i,1)-data(i,2));
                min_index=i;
            end
        end
        
    end
    
    if(delta>abs(data(i+1,1)-data(i+1,2)))
        o_data=[o_data;data(i+1,:)];
    else
        o_data=[o_data;data(min_index,:)];
    end
    
    res=o_data;
end



