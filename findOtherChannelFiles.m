function dirIndices = findOtherChannelFiles(D, FileName)
    %D = dir
    dirIndices=0;
    strFileName=strsplit(FileName,'_');
    strFileName=strcat(strFileName(2:end));    %get string after ch #
            for i=1:length(D)
                dd=strsplit(D(i).name,'_'); %get string after ch #
                if (length(dd)-1)==length(strFileName)
                    cmpFileNames=strcmp(dd(2:end),strFileName);
                    sameFile=prod(cmpFileNames);              %multiply values of matrix, if 1 then rest of file name is same and is correspondign channel for same data
                    %if isempty(cmpFileNames)
                    %    sameFile=0;
                    %end
                    % d=strfind(D(i).name,FileNum);     %find other channel files with same FileNum, ie. '008'
                    if sameFile       %when filename contains same string after Channel # this is another channel matching to same measurement
                        if dirIndices==0
                            dirIndices=i;       %store first index of channel data num
                        else
                            dirIndices=[dirIndices i];
                        end

                    end
                end
            end
end