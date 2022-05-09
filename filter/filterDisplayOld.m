function filtered = filterDisplay(samplesperscan, fc, data)
filtered=data; 
samplesperscan=size(data,1);

switch samplesperscan
    case 175
        filtered = imgaussfilt(data,[2.5 2]);             %fo r 4 channel %changed from 2 to 3.5 08/14/2017
    case 500
        filtered = imgaussfilt(data,[8 2]);         %for 2 channel %changed from 10 08/16/2017
    case 1000
        filtered = imgaussfilt(data,[16 2]);  %changed from 20 08/16/2017 for 1 channel
    otherwise
        filtered = imgaussfilt(data,[3 2]);       %default%changed 08/14/2017 %changed from 4 08/16/2017
end

%filtered=data; 
end