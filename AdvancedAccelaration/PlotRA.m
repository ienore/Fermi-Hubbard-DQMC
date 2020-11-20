function data_plot = PlotRA(temp_observable,RA)
    back_track = round(length(temp_observable)/RA);
    data_plot = zeros([1,length(temp_observable)-back_track]);
    for index_b =  1:1:length(temp_observable)-back_track
        if(index_b == 1)
            temp_data = temp_observable(index_b:index_b+back_track);
            data_plot(index_b) = mean(temp_data);
            s = sum(temp_data);
        else
            s = s - temp_observable(index_b-1)+temp_observable(index_b+back_track);
            data_plot(index_b)=s/back_track;
        end
    end
end

