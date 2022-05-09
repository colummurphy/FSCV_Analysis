function [Isub,BGavgi,BGrefmatrix]=refreshPlots(Iread,sampling_rate,BGstartpoint,BGavg,t_end)
    %refreshPlots
    %initialize original plot based on loaded data
    Iread_total=Iread;          %for averaging later as option
    Imeas=Iread;
    %Imeas=Iread';
    %Imeas=flipud(Imeas);    %flip so goes from negative to positive to negative potential
    sizeMeas=size(Imeas);

    %set time scale
    t_range=sizeMeas(2);        %total # timepoints
    orig_t_range=t_range; %display(orig_t_range)
    if t_end>t_range
        t_end=t_range;
    end

    %plot current subtracted from avg of ##BGavg points from BGstartpoint
    if BGstartpoint<=0
        BGstartpoint=BGavg+1;
    end
    if (BGstartpoint-BGavg)<=0
        BGstartpoint=BGstartpoint+BGavg;
    end
    BGend=(round(BGstartpoint)+BGavg);
    if BGend>t_range
        BGend=t_range;
    end
    BGvectors=Imeas(:,(round(BGstartpoint)-BGavg):BGend); 
    BGavgi=mean(BGvectors,2);         %mean ## sample points BG for subtraction
    BGrefmatrix=repmat(BGavgi,1,orig_t_range);   %tile BGavg vector for entire meas matrix span
    Isub=Imeas-BGrefmatrix;     %subtracted current data
    %PCA components initialize to subtracted current data
    EPH_display=Isub;       %initialize default;
    E_display=Isub;
    EDA_display=Isub;
    
end