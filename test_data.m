input_data_ns=load('data/9_6/p_sweep_grating_250_5_350_n=3.5_d=200_1:10_nostruct.dat');
input_data_ws=load('data/9_6/p_sweep_grating_250_5_350_n=3.5_d=200_1:10.dat');
input_data_ws(:,4)=input_data_ws(:,4)./input_data_ns(:,4);

t=0;
all_match_data=[];
%input_data=load('data/p_sweep_grating_250_5_350_n=3.5_d=200_1:10.dat');
input_data=input_data_ws;
data_size=size(input_data);
col_length=data_size(2);
[out_data, cnt]=split_data(input_data,1);
[dispersion_data,c_t]=numeric_dispersion([3.5 1 1], 200e-6);
t=t+c_t;


% figure
for i=1:cnt
    data=out_data(:,col_length*(i-1)+1:col_length*i);
    grating_period=data(1,1);
    c_t=draw_curve(dispersion_data,data,grating_period,[3.5 1],'9_6');
    t=t+c_t;
    [ p_omiga, beta,m,n, c_t ] = match_cal_grating_beta( grating_period, dispersion_data,2,'9_6');
    t=t+c_t;
    [ match_data,c_t] = match_cal_sim_data(data, p_omiga,m,n);
    %disp(match_data);
    t=t+c_t;
    if(~isempty(match_data))
        h=figure(1);
        plot(data(:,2),data(:,3));
        hold on
        plot(match_data(:,2)*1e12,match_data(:,3),'*r');
        plot(data(:,2),data(:,4)/max(data(:,4)),'--k');    
        for m=1:3
            g_omiga=2*pi*3e8*m/grating_period/3.5/2/pi;
            plot(ones(1,length(data))*g_omiga,linspace(0,1,length(data)),'-b');
        end
        hold off
        len=size(match_data);
        len=len(1);
        for j=1:len
              
            str = sprintf('(%d,%d)',match_data(j,5),match_data(j,6));
            text(match_data(j,2)*1e12,match_data(j,3),str);
            
        end
        str_path=sprintf('pic/9_6/p=%e.png',grating_period);
        saveas(h,str_path);
        close(h)
    %     if ~isempty(match_data)
        all_match_data = [all_match_data; [ones(length(match_data(:,2)),1)*grating_period,match_data(:,2)]];
    %     end
    end
%     plot(data(:,2),data(:,3));
%     hold on
end
% hold off
figure
plot(all_match_data(:,1),all_match_data(:,2),'.');
