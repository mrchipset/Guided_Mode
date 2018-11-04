function [ output_data,cnt ] = split_data( input_data,index )
%split_data Summary of this function goes here
%   Detailed explanation goes here
    tmp_data=input_data(:,index);
    diff_data=diff(tmp_data);
    flag_index=find(sign(diff_data));
    cnt=length(flag_index)+1;
    output_data=[];
    output_data=[output_data;input_data(1:flag_index(1),:)];
    for i=1:length(flag_index)
        

        if(i==length(flag_index))
            output_data=[output_data input_data(flag_index(i)+1:end,:)];
        else
            output_data=[output_data input_data(flag_index(i)+1:flag_index(i+1),:)];
        end
    end
    %output_data=[output_data;input_data(flag_index(i):end,:)];
        
end

