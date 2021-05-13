clc;
clear;
close all;
global TURN_ON_PLOTS;
TURN_ON_PLOTS = true; % set this to true if you want the figures to appear.
global REMOVE_CHANNEL_SIX;
% Surprisingly, we get consistently worse performance when we remove channel 6. There's some good information there!
REMOVE_CHANNEL_SIX = false; 

fdsGroup4 = fileDatastore(fullfile('LDACLASSIFYG4'), 'PreviewFcn', @load, 'ReadFcn', @load, 'IncludeSubfolders', false, 'FileExtensions', '.mat');
previewData = preview(fdsGroup4); % Peeks into the first file in the data directory
fs=previewData.samplingRate; %sampling rate / frequency. (~3004 Hz)
dt=1/fs;

numberOfFiles = length(fdsGroup4.Files);
% Load the files, pulling their names from the filename.
for idx = 1:numberOfFiles
    [~, name, ~] = fileparts(fdsGroup4.Files{idx});
    mymodel.name{idx} = name;
    emgDataStruct = read(fdsGroup4);
    mymodel.data{idx} = emgDataStruct.Data;
    if REMOVE_CHANNEL_SIX
        mymodel.data{idx}(6) = []; %channel 6 might be bad. Remove?
    end
    mymodel.length{idx} = emgDataStruct.length_sec;
end
numberOfChans= length(mymodel.data{1});%%find how many channels we have

%% Construct the filters, and run data through them.
% We used a single bandpass filter with the filtfilt function, as well as
% spectrum interpolation to remove the 60 Hz + harmonic noise. Please note
% that it is impossible to get a traditional filter plot from the 60 Hz &
% harmonic filtering we did! Hence why we only show the before and after
% plots after combining both the bandpass and spectrum interpolation.

lp=480; % Hz
hp=30; % Hz

bp1  = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp, ...
    'SampleRate',fs);

fvtool(bp1);
title("Bandpass filter, passband = 30-480 Hz");

total_data = cell(1,numberOfFiles);
for kk=1:numberOfFiles
    data=mymodel.data{kk};
    data=cell2mat(data); % converts data into matrix format.

    % band pass
    fulldata = filtfilt(bp1, data); % No phase shift.
    % Spectrum interpolation to remove 60 hz + harmonic noise. Works in the
    % frequency domain and does not rely on filters. No phase shift.
    [h,w]=size(data);
    for idx = 1:w
        fulldata(:, idx) = spectrumInterpolation(fulldata(:, idx), fs, 60, 3, 2);
    end
    
    % plot code
    if TURN_ON_PLOTS
        for channels=1:width(data)
            f = figure;
            nexttile;
            plot(data(:,channels));
            title(['Raw EMG, Ch ', num2str(channels)]);
            nexttile;
            plot(fulldata(:,channels));
            title(['Filtered EMG, Ch ', num2str(channels)]);
            nexttile;
            pwelch(data(:,channels), [], [], [], fs);
            title(['Raw EMG PSD Ch ', num2str(channels)]);
            nexttile;
            pwelch(fulldata(:,channels), [], [], [], fs);
            title(['Filtered EMG PSD Ch ', num2str(channels)]);
            title(f.Children, convertCharsToStrings(mymodel.name{kk}));
        end
    end
    
    mymodel.data{kk}=fulldata; % write back to conserve memory.
    mymodel.meanChan{kk} = mean(fulldata,2);
end

save('filteredData.mat','mymodel','numberOfFiles','numberOfChans','fs','dt')

%% This section of code removes the off data and places the "on" data into the onData member within the mymodel structure.
% There is also a member called "overlayData" that one can use to make nice
% plots to verify that the onData seems correct.

mymodel = removeOffData(mymodel);
% We need to make the datasets all the same length now.
% Determine the shortest onData length, then just shorten all the other
% datasets to match it.

shortestPose = length(mymodel.onData{1});
for i=1:numberOfFiles
    if(shortestPose > length(mymodel.onData{i}))
        shortestPose = length(mymodel.onData{i});
    end
end

for i=1:numberOfFiles
    tempDat = mymodel.onData{i};
    tempDat = tempDat(1:shortestPose, :);
    mymodel.onData{i} = tempDat;
end

TrimmedTF = mymodel.onData;

%% Feature extraction
% Here, we determine the binning and 4 features (MAV, SSC, WL, and ZC) used
% to create our LDA classifier.
binTime = 100; % milliseconds
binsize=(binTime/1000)*fs;

Tz = cell(numberOfFiles, 1);
for kk=1:numberOfFiles
    Tz{kk} = TrimmedTF{kk}'; % Feature extraction requires our data to be transposed, so we'll do that here.
end

% extract each of our 4 choosen features, please see respective files
for jj=1:numberOfFiles
    nBin=floor(length(Tz{jj})/binsize); %%how many bins do we have 
    Bz=floor(linspace(1,length(Tz{jj}),nBin)); %%the bins themselves
    Ez=[];%%an empty matrix for intiallizing the other
    ZC_TZ=Ez; %%zero crossings
    SSC_TZ=Ez;%%slope sign changes
    MAV_TZ=Ez;%%mean absolute value
    WL_TZ=Ez;%% wavlength
    
    for kk=1:numberOfChans    
        for ii=1:nBin-1
            Dz=Bz(ii+1);%%what bin we are in
            
            ZC_TZ(kk,ii)=ZCz(Tz{jj}(kk,Bz(ii):Dz));%%find zero crossing
            MAV_TZ(kk,ii)=MAVz(Tz{jj}(kk,Bz(ii):Dz));%%find mean absolute value
            SSC_TZ(kk,ii)=SSCz(Tz{jj}(kk,Bz(ii):Dz));%%find slope sign change
            WL_TZ(kk,ii)=WLz(Tz{jj}(kk,Bz(ii):Dz)); %%find wavelngth
        end
    end
    
    ZC{jj}=ZC_TZ;%%add the current ZC feature to the structure
    MAV{jj}=MAV_TZ;%%add the current MAV feature to the structure
    SSC{jj}=SSC_TZ;%%add the current SSC feature to the structure
    WL{jj}=WL_TZ;%%add the current WL feature to the structure
end

% Plot the features
if TURN_ON_PLOTS
    for idx=1:numberOfFiles
        for jdx=1:numberOfChans
            f = figure;
            nexttile;
            plot(ZC{idx}(jdx, :));
            title(['Zero Crossing, Ch ', num2str(jdx)]);
            nexttile;
            plot(MAV{idx}(jdx, :));
            title(['Mean Absolute Value, Ch ', num2str(jdx)]);
            nexttile;
            plot(SSC{idx}(jdx, :));
            title(['Slope Sign Change, Ch ', num2str(jdx)]);
            nexttile;
            plot(WL{idx}(jdx, :));
            title(['Waveform-Length, Ch ', num2str(jdx)]);
            title(f.Children, [mymodel.name{idx}, ' Features. ', num2str(nBin-1), ' bins, length = ', num2str(binTime), ' ms.']);
        end
    end
end

%Data separation
trainingIndex = 1:2:24;
validationIndex = 2:2:24;

for i = 1:12
    indexT = trainingIndex(i);
    trainingNames{i} = mymodel.name{indexT};

    indexV = validationIndex(i);
    validationNames{i} = mymodel.name{indexV};
end

ZC_train = ZC(1, trainingIndex);
MAV_train = MAV(1, trainingIndex);
SSC_train = SSC(1, trainingIndex);
WL_train = WL(1, trainingIndex);

ZC_val = ZC(1, validationIndex);
MAV_val = MAV(1, validationIndex);
SSC_val = SSC(1, validationIndex);
WL_val = WL(1, validationIndex);

numberOfPoses = numberOfFiles / 2;

%% Matrices used for LDA
% Create the class mean matrix, and the associated between class separation
% matrix, as well as the sorted matrix of eigenvectors (Wsroted), thresheld above
% 1e-6. This matrix will serve as our linear transformation from the
% original feature space to the newly defined space.

classMeansMatrix = zeros(numberOfChans * 4, numberOfPoses);
for i=1:numberOfPoses
    featureMatrix =[WL_train{i};SSC_train{i};MAV_train{i};ZC_train{i}];%create a matrix of our feature
    % For each posture set, we will determine the average value of each
    % feature extracted across the number of sample bins.
    classMeansMatrix(:,i) = mean(featureMatrix,2)';
end
betweenClassMatrix = cov(classMeansMatrix'); % again, we need to take the transpose to get the proper 32x32 matrix out from the covariance calculation.

allClassMeanVector = mean(classMeansMatrix); % This is just the mean of all class vectors. Don't really need this, but it's nice to have just in case.

% Create the within class scatter matrix, which is just the average of each
% output classes' covariance matrices.

withinClassMatrix = zeros(numberOfChans * 4,numberOfChans * 4);
for i=1:numberOfPoses
    featureMatrix =[WL_train{i};SSC_train{i};MAV_train{i};ZC_train{i}];%create a matrix of our feature
    featureMatrix = featureMatrix'; % Transpose it to properly compute the covariance matrix for this class.
    withinClassMatrix =  withinClassMatrix + cov(featureMatrix);
end

withinClassMatrix = withinClassMatrix / numberOfPoses;


% Now that we have the prereq matrices, we can create the optimizing
% criterion matrix, which is the crux of the LDA algorithm.

optMatrix = withinClassMatrix \ betweenClassMatrix; % This is equivalent to inv(A) * b
[W,Z]=eig(optMatrix); % W is a matrix of eigenvectors, Z is a matrix of eigenvalues.

% We should now sort the eigenvalues and their associated eigenvectors.
[z, ind] = sort(diag(Z), 'descend');
Zsorted = Z(ind, ind);
Wsorted = W(:, ind);

% Now we can remove eigenvalues and the eigenvectors that are very small
% (<1e-6)
Wsorted = Wsorted(:, 1:11); % After the 10th eigenvalue, these vectors get really small.


%% Euclidean Distance
% This calculation determines each data points' distance from the mean
% vectors in our newly defined space. The mean vector that the point is
% closest to is the output of the LDA classifier.

% Now that we have our matrix, we can go ahead and transform the validation
% data into our newly defined space, and we also need to transform the
% class mean vectors.
transformedClassMeanMatrix = real(Wsorted)' * classMeansMatrix;

for i=1:numberOfPoses
    a=[WL_val{i};SSC_val{i};MAV_val{i};ZC_val{i}];
    Yval{i} = real(Wsorted)' * a;
    Y_avg{i} = transformedClassMeanMatrix(:,i);
end

% calculate euclidean distance from each point in validation data to class
% average values for each training class
% result is structure with 12 cells of 12xnumBin matrices
% each cell corresponds to one validation data pose
% each row in the cells represents distance to a training data pose
for m = 1:numberOfPoses         % val poses
    for n = 1:numberOfPoses     % train poses
        numBins = length(Yval{n});
        clear dist sqTerms
        
        for i = 1:numBins
            for j = 1:width(Wsorted)
                yVal_i = Yval{1,m}(j,i);
                yTrain_i = Y_avg{1,n}(j);
                sqTerms(j) = (yVal_i - yTrain_i)^2;
            end            
            dist(i) = sqrt(sum(sqTerms));
        end
        distByTrainPose(n,:) = dist;
    end
    distStruct{m} = distByTrainPose; 
end

%% Classification
% This is where the confusion matrix is created and we can view the results
% of our classification.
classNames = ['Grasping' 'Hand Close' 'Hand Open' 'Off' 'Thumb Abduction'...
    'Thumb Adduction' 'Wrist Extension' 'Wrist Flexion' 'Wrist Pronation'...
    'Wrist Rad Dev' 'Wrist Supination' 'Wrist Ulnar Dev'];
outPutFull=[];
for i = 1:numberOfPoses
    valPose = distStruct{i};
    min_d=10^6;%%we need to set a treshold
    for j=1:numBins
       valBin=valPose(:,j);
       [minDist,index] = min(valBin);
       outputClass(j) = index;
    end
    outPutFull=[outPutFull outputClass];
end
outPutFull = categorical(outPutFull);
trueCats=[];
for ii=1:numberOfPoses
    trueCats=[trueCats ones(1,nBin-1)*ii];
end
trueCats=categorical(trueCats);
plotconfusion(trueCats,outPutFull)
if REMOVE_CHANNEL_SIX
    % We actually get worse performance with channel 6 removed.
    title(['Confusion Matrix with Channel 6 removed, binsize = ', num2str(binTime), ' milliseconds']); 
else
    title(['Confusion Matrix with Channel 6 intact, binsize = ', num2str(binTime), ' milliseconds']);
end
