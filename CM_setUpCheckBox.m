function [] = CM_setUpCheckBox(checkBoxes, checkNums)
global plotParam

%{
  Function accepts an array of check boxes and an array
  of check numbers.

  Sets the check box label text and the label text color. 
  Sets the check box initial value.
%}

if length(checkBoxes) ~= length(checkNums)
    return;
end    

numBoxes = length(checkBoxes);
for i = 1:numBoxes

    % get check box and number
    checkBox = checkBoxes(i);
    checkNum = checkNums(i);
      
    % set color
    checkBox.FontColor = plotParam.colorFSCV(checkNum,:);
            
    % set label    
    checkLabel = ['ch' num2str(checkNum)];           
    if ~isempty(plotParam.sites)
        checkLabel = plotParam.sites{checkNum};
    end
    checkBox.Text = checkLabel;

    % if checkNum is in selch array, plot the channel
    plotValue = 0;
    if ismember(checkNum, plotParam.selch)
        plotValue = 1;
    end  
    checkBox.Value = plotValue;

end

end