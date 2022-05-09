function filtered = filterDisplay(samplesperscan, fc, data)
%changed 06/23, because imgaussfilt not available on matlab 2013 on chunky
filtered=data; 
samplesperscan=size(data,1);

switch samplesperscan
    case 175
        %ff = fspecial('gaussian', [2.5 2], 1);
        if exist('imgaussfilt')>0
                %if is function, not function  in Matlab 2013 and earlier
                %versions                
            filtered = imgaussfilt(data,[2.5 2]);             %fo r 4 channel %changed from 2 to 3.5 08/14/2017

        else
            sigy=2.5;
            sigx=1;
            %sizey=2*ceil(2*sigy)+1;
            %sizex=2*ceil(2*sigx)+1;
            sizey=20;
            sizex=7;
            sig=(sizey-1)/2;
            ff = fspecial('gaussian', [sizey sizex], 2);    
            filtered = filter2(ff, data);
        end
    case 500
        filtered = imgaussfilt(data,[8 2]);         %for 2 channel %changed from 10 08/16/2017
    case 1000
        filtered = imgaussfilt(data,[16 2]);  %changed from 20 08/16/2017 for 1 channel
    otherwise
        filtered = imgaussfilt(data,[3 2]);       %default%changed 08/14/2017 %changed from 4 08/16/2017, for rodent, 214 samples
end

%filtered=data; 
end