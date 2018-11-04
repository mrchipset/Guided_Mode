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
    for j = i.beta
        omiga=[omiga i.omiga];
        beta=[beta j];
        
    end

    p=[p 2*pi/300e-6];
    p_2=[p_2 4*pi/300e-6];
    p_3=[p_3 6*pi/300e-6];
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

for i=[2*pi/300e-6 4*pi/300e-6 6*pi/300e-6]
    d_beta=beta((abs(beta-i)/1e4)<0.01);
    d_omiga=omiga((abs(beta-i)/1e4)<0.01);
    diff(d_omiga);
    d_beta=d_beta(diff(d_omiga)>0.1*1e12);
    d_omiga=d_omiga(diff(d_omiga)>0.1*1e12);
    p_beta=[p_beta d_beta];
    p_omiga=[p_omiga d_omiga];
    
%     d_omiga_flag=diff(d_omiga)>0.1*1e12;
%     d_omiga_flag=diff(d_omiga)>0.1*1e12;
%     d_beta_flag=d_beta(diff(d_omiga)>0.1*1e12);
%     d_omiga_flag=diff(d_omiga)>0.1*1e12;
%     index=0;
%     cnt=0;
%     sum_omiga=0;
%     sum_beta=0;
%     for j=d_omiga_flag
%         index=index+1;
%         cnt=cnt+1;
%         sum_omiga=sum_omiga+d_omiga(index);
%         sum_beta=sum_beta+d_beta(index);
%         if(j==1)||index==length(d_omiga_flag)
%             p_omiga=[p_omiga sum_omiga/cnt];
%             p_beta=[p_beta sum_beta/cnt];
%             cnt=0;
%             sum_omiga=0;
%             sum_beta=0;
%         end
%     end
        
end

%p_beta=sort(p_beta);
%p_omiga=sort(p_omiga);

% delete(pl);
omiga=omiga/2/pi;
p_omiga=p_omiga/2/pi;
c=c/2/pi;

plot(beta,omiga,'.k');
%plot(omiga,beta,'.k');
hold on

p=ones(100,1);
p_2=ones(100,1);
p_3=ones(100,1);
p=p*2*pi/300e-6;
p_2=p_2*4*pi/300e-6;
p_3=p_3*6*pi/300e-6;
plot(p,linspace(0,max(omiga),length(p)));
plot(p_2,linspace(0,max(omiga),length(p_2)));
plot(p_3,linspace(0,max(omiga),length(p_3)));
plot(beta,beta*c/n1)
plot(beta,beta*c/n2)
plot(p_beta,p_omiga,'*r');

hold off
axis([0 7e4 0 1e12])