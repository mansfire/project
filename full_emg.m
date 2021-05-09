clc;
clear;
close all;

fdsGroup4 = fileDatastore(fullfile('LDACLASSIFYG4'), 'PreviewFcn', @load, 'ReadFcn', @load, 'IncludeSubfolders', false, 'FileExtensions', '.mat');
previewData = preview(fdsGroup4); % Peeks into the first file in the data directory
fs=previewData.samplingRate ;%sampling rate
dt=1/fs;

numberOfFiles = length(fdsGroup4.Files);
% Load the files, pulling their names from the filename.
for idx = 1:numberOfFiles
    [~, name, ~] = fileparts(fdsGroup4.Files{idx});
    mymodel.name{idx} = name;
    emgDataStruct = read(fdsGroup4);
    mymodel.data{idx} = emgDataStruct.Data;
    %mymodel.data{idx}{6} = zeros(length(mymodel.data{idx}{6}), 1); %channel 6 might be bad. Remove?
    mymodel.length{idx} = emgDataStruct.length_sec;
end
numberofChans= length(mymodel.data{1});%%find how many channels we have
%% Construct the filters, and run data through them.
lp=480; % Hz
hp=30; % Hz

bp1  = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp, ...
    'SampleRate',fs);

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
    mymodel.data{kk}=fulldata; % write back to conserve memory.
    mymodel.meanChan{kk} = mean(fulldata,2);
end

%% add all the channels for each of the 11 test and validation postures
for ii=1:numberOfFiles
    total_data{idx}=sum(abs(mymodel.data{idx}), 2); % Rectifies the data, then adds all channels together.
end


%% This section of code removes the off data and places the "on" data into the onData member within the mymodel structure.
% There is also a member called "overlayData" that one can use to make nice
% plots to verify that the onData seems correct.
mymodel = removeOffData(mymodel);
TrimmedTF = mymodel.onData;

% filename='trimmed.mat';
% save(filename)
%% off data prep
binsize=0.05*fs;

Tz = cell(numberOfFiles, 1);
for kk=1:numberOfFiles
    Tz{kk} = TrimmedTF{kk}'; % Feature extraction requires our data to be transposed, so we'll do that here.
end

% extract each of our 4 choosen features, please see respective files
for jj=1:numberOfFiles
    nBin=floor(length(Tz{jj})/binsize); %%how many bins so we have 
    Bz=floor(linspace(1,length(Tz{jj}),nBin)); %%the bins themselves
    Ez=[];%%an empty matrix for intiallizing the other
    ZC_TZ=Ez; %%zerocrossings
    SSC_TZ=Ez;%%slope sign changes
    MAV_TZ=Ez;%%mean absolute value
    WL_TZ=Ez;%% wavlength
    for kk=1:numberofChans
        
        for ii=1:nBin-1
            Dz=Bz(ii+1);%%what bin we are in
            
            ZC_TZ(kk,ii)=ZCz(Tz{jj}(kk,Bz(ii):Dz));%%find zero corrsing
            MAV_TZ(kk,ii)=MAVz(Tz{jj}(kk,Bz(ii):Dz));%%find mean absolute value
            SSC_TZ(kk,ii)=SSCz(Tz{jj}(kk,Bz(ii):Dz));%%find slope sign change
            WL_TZ(kk,ii)=WLz(Tz{jj}(kk,Bz(ii):Dz)); %%find wavelngth 
            
        end
        ZC{jj}=ZC_TZ;%%add the current ZC feature to the structure
        MAV{jj}=MAV_TZ;%%add the current MAV feature to the structure
        SSc{jj}=SSC_TZ;%%add the current SSC feature to the structure
        WL{jj}=WL_TZ;%%add the current WL feature to the structure
        
    end
    
end
cord=0;
%% covariance
for ii=1:numberOfFiles
    a=[WL{ii};SSc{ii};MAV{ii};ZC{ii}];%create a matrix of our feature
    cord=cord+1; %%incriminent the index
    [M,N]=size(a);
    mn=mean(a,2);
    a=a-repmat(mn,1,N);
    cov_mat{cord}=cov(a*a');%%create the covariance natrix
    a_mat{ii}=a;%%save feature

end
for ii=1:numberOfFiles
    [V,D] = eig(cov_mat{ii});
    D=diag(D);%%look at the acual values
    [val,index]=sort(-1*D);
    index=index(1:3);%%only keep the best

    D=D(index);
    V=V(:,index);
    eig_vec{ii}=V;%%eigne vector
end

for ii=1:2:numberOfFiles
        euc=0;
        min_distance=10^12;
        for jj=1:length(index)
        meandist=mean(W{ii}(jj,:));
        euc=euc+dist(meandist,W{ii+1}(jj,:));


end
end
