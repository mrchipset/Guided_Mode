function [ omiga, beta, r_m, r_n, c_t ] = match_cal_grating_beta( grating_period, data, level, path_prefix)
%match_cal_grating_beta Summary of this function goes here
%   Detailed explanation goes here
    %%
    %Dump the data
    tic;
    omiga=[];
    beta=[];
    d_beta=[];
    d_omiga=[];
    
    for i=data  
            omiga=[omiga i.omiga];
            beta=[beta i.beta];
    end
    
    p_beta=[];
    p_omiga=[];
    p_data=[];
    f=fopen('log/grating_cal.txt','w');
    fprintf(f,'Grating Beta,\tCalculation Beta,\tCalculation Freq\n');
    for j=2:2:2*level;
        i=j*pi/grating_period;
        d_beta=beta((abs(beta-i)/i)<0.01);
        d_omiga=omiga((abs(beta-i)/i)<0.01);
    %     d_beta=d_beta(diff(d_omiga)>0.1*1e12);
    %     d_omiga=d_omiga(diff(d_omiga)>0.1*1e12);
        p_beta=[p_beta d_beta];
        p_omiga=[p_omiga d_omiga];
        
        for l=1:length(d_beta)
            fprintf(f,'%.4f,\t%.4f,\t%.4f\n',i,d_beta(l),d_omiga(l));
            p_data=[p_data;i,d_beta(l),d_omiga(l)];
        end
        

    end
    fclose(f);

    p_data=drop_data1(p_data);
%     p_omiga=p_data(:,3);
%     p_beta=p_data(:,1);
    p_omiga=[];
    p_beta=[];
    p_m=[];
    p_n=[];
    str_path=sprintf('log/%s/data_dropped_grating_cal_p=%e.txt',path_prefix, grating_period);
    %f=fopen('log/data_dropped_grating_cal.txt','w');
    f=fopen(str_path,'w');
    fprintf(f,'Grating Beta,\tCalculation Beta,\tCalculation Freq,\tm,\tn\n');
    m=0;
    n=0;
    tmp=0;
    for i=1:length(p_data)    
        if(tmp~=p_data(i,1))
            m=m+1;
            n=0;
            tmp=p_data(i,1);
        end
        n=n+1;
        if(n<=3)
            p_omiga=[p_omiga;p_data(i,3)];
            p_beta=[p_beta;p_data(i,1)];
            p_m=[p_m;m];
            p_n=[p_n;n];
            fprintf(f,'%.4f,\t%.4f,\t%.4f,\t%d,\t%d\n',p_data(i,1),p_data(i,2),p_data(i,3)/1e12/2/pi,m,n);
        end
    end
    fclose(f);
    omiga=p_omiga;
    beta=p_beta;
    r_m=p_m;
    r_n=p_n;
    c_t=toc;

end


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


