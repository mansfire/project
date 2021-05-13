function newStruct = removeOffData(inputStructure)
% These values were handcrafted, and are tailored to each posture / trial.
% To reproduce, just use rmData(inputStructure.data{idx}, cutOff, true,
% channel)
global REMOVE_CHANNEL_SIX;
newStruct = inputStructure;
channelSelection = [4 2 3 5 8 8 1 1 8 1 7 8 8 8 4 4 6 6 1 1 2 4 8 8];

if REMOVE_CHANNEL_SIX
    for i=1:length(channelSelection)
        if channelSelection(i) >= 6
            % just move everything down by a channel.
            channelSelection(i) = channelSelection(i) - 1;
        end
    end
end
cutoffSelection = [0.3 -0.2 -0.3 -0.3 0 -0.2 -100 -100 0.6 0.2 0.2 0.1 0.1 0.1 -0.3 0 -0.4 0 -0.4 -0.4 -0.2 -0.2 -0.2 -0.2];

for i=1:length(newStruct.data)
    [newStruct.onData{i}, newStruct.overlayData{i}] = rmData(newStruct.data{i}, newStruct.name{i}, cutoffSelection(i), false, channelSelection(i));
end
end


function [scrubbedData, overlayData] = rmData(inputData, poseName, threshold, turnOnDebugPlots, channelToUse)
% This function is intended to remove off data. It will use the threshold
% std dev to determine the data that would be considered "on".
% Input parameters:
% inputData - matrix of the data, with each column representing a channel.
% threshold - multiple of the standard deviation
% turnOnDebugPlots - Used to determine the correct value for the threshold.
% channelToUse - The channel to use for data selection. Ignored when
% turnOnDebugPlots is True.
% Output - scrubbedData - matrix of the data, with the same number of
% columns as inputData, but with off sections removed.
global TURN_ON_PLOTS;
smoothData=smoothdata(abs(inputData),'sgolay',round(length(inputData) / 100));
dataMean=mean(smoothData);
dataStd=std(smoothData);
calculatedThreshold = dataMean + (threshold * dataStd);

if turnOnDebugPlots == true
    figure;
    t = 1:length(inputData);
    for i=1:width(smoothData)
        nexttile;
        plot(smoothData(:,i))
        hold on;
        y = ones(length(t),1) * calculatedThreshold(i);
        plot(t, y)
        hold off;
        title(['Channel ', num2str(int32(i)), ', threshold = \mu + ', num2str(threshold), '\sigma ']);
    end
else
    % If a value passes above the calculated threshold, it should be considered
    % as "on" data, we'll use the best looking channel as determined by the
    % debug plot above.
    
    % Create a logical array with the same length as the data, then use it as
    % a mask for the rest of the data.
    
    channelData=smoothdata(abs(inputData(:, channelToUse)),'sgolay',round(length(inputData) / 100));
    channelData(channelData<=calculatedThreshold(channelToUse)) = 0;
    channelData(channelData>calculatedThreshold(channelToUse)) = 1;
    maskArray = logical(channelData);
    
    % It is important to not use smoothData here as the next step of this
    % process is to do feature extraction. Smooth data was only used to get
    % the "waveform" of the emg data for easier off data removal.
    
    % Creating this just to get the size of scrubbedData.
    chanData = inputData(:,1);
    chanData = chanData(maskArray);
    scrubbedData = zeros(length(chanData), width(inputData));
    overlayData = zeros(size(inputData));
    for i=1:width(inputData)
        chanData = inputData(:,i);
        scrubbedData(:,i) = chanData(maskArray);
        overlayData(:,i) = chanData .* channelData;
    end
    % Uncomment this section if you want to verify that the "on data" looks
    % right. This will create a lot of plots!!! Be warned.
    if TURN_ON_PLOTS
        for i=1:width(inputData)
            f = figure;
            nexttile;
            plot(inputData(:,i))
            hold on;
            plot(overlayData(:,i));
            hold off;
            title(['On Data Overlaid on Filtered Signal']);
            nexttile;
            plot(scrubbedData(:,i));
            title(['On Data Only']);
            title(f.Children, [poseName, ' Off Data Removal Plots, Ch ', num2str(i)']);
        end
    end
end

end