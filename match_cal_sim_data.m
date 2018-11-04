function [ match_data, t ] = match_cal_sim_data( sim_data, p_omiga, m, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    tic;
    p_omiga=p_omiga/2/pi;
  
    
    S11=sim_data;

    freq=S11(:,2);
    ref=S11(:,3);
    intens=S11(:,4);
    %intensity=S11(:,4)/max(S11(:,4));
%     tar_ref=[];
%     tar_freq=[];
%     tar_intens=[];
%     for i=-1:1:1
%         tar_ref=[tar_ref ref(find(diff(sign(diff(ref)))~=0)-i)];
%         tar_freq=[tar_freq freq(find(diff(sign(diff(ref)))~=0)-i)];
%         tar_intens=[tar_intens intens(find(diff(sign(diff(ref)))~=0)-i)];
%     end
    tar_ref=ref(diff(sign(diff(ref)))~=0);
    tar_freq=freq(diff(sign(diff(ref)))~=0);
    tar_intens=intens(diff(sign(diff(ref)))~=0);
    len=length(tar_freq);
    p_omiga=p_omiga/1e12;
    tar_freq=tar_freq/1e12;


    res=[];
    f=fopen('log/cal_sim_point.txt','w');
    fprintf(f,'Calculation Freq,\tSimulation Freq,\tSimulation Reflection\n');
    c_data=[];
    for id=1:length(p_omiga)
        i=p_omiga(id);
        for index=1:len
            if(abs(tar_freq(index)-i)/i<0.05||abs(tar_freq(index)-i)/tar_freq(index)<0.05)
    %             fprintf(f,'%.4f,\t%.4f,\t%.4f\n',i,tar_freq(index),tar_ref(index));
    %             plot(tar_freq(index)*1e12,tar_ref(index),'*r');
                c_data=[c_data;i,tar_freq(index),tar_ref(index),tar_intens(index),m(id),n(id)];
            end
        end
    end

    c_data=drop_data2(c_data);
    len=size(c_data);
    len=len(1);
    for i = 1:len
        fprintf(f,'%.4f,\t%.4f,\t%.4f\n',c_data(i,1),c_data(i,2),c_data(i,3));
    end
    fclose(f);
    match_data=match_intensity(c_data,sim_data,1.3);
    t=toc;


end
function r_data=match_intensity(o_data,data,intens_level)

    d_data=[];
    tmp_o_data=o_data;
    tmp=data;
    max_p=max(tmp(:,4));
    tmp(:,4)=tmp(:,4)/max(tmp(:,4));
    tmp_o_data(:,4)=tmp_o_data(:,4)/max_p;
    intensity=tmp(:,4);
    intensity=sort(intensity);
    max_int=intensity(end-3:end);
    len=size(tmp_o_data);
    len=len(1);
    for i=1:len
        for j=length(max_int)            
            if(abs(tmp_o_data(i,4))>mean(intensity)*intens_level)
                d_data=[d_data;o_data(i,:)];
            end
        end
    end
    
    r_data=d_data;
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