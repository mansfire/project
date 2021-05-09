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
    for idx = 1:width(data)
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

%% off data prep
binsize=0.05*fs;

for kk=1:numberOfFiles
    clear T;
    for ii=1:8
        
        for jj=1:length(TrimmedTF{kk}(:,ii))
            T(ii,jj)=(TrimmedTF{kk}(jj,ii));
            
        end
    end
    Tz{kk}=T;
end
% extract each of our 4 choosen features, please see respective files
for jj=1:numberOfFiles
    nBin=floor(length(Tz{jj})/binsize);
    Bz=floor(linspace(1,length(Tz{jj}),nBin));
    Ez=[];
    ZC_TZ=Ez;
    SSC_TZ=Ez;
    MAV_TZ=Ez;
    WL_TZ=Ez;
    for kk=1:8
        
        for ii=1:nBin-1
            Dz=Bz(ii+1);
            
            ZC_TZ(kk,ii)=ZCz(Tz{jj}(kk,Bz(ii):Dz));
            MAV_TZ(kk,ii)=MAVz(Tz{jj}(kk,Bz(ii):Dz));
            SSC_TZ(kk,ii)=SSCz(Tz{jj}(kk,Bz(ii):Dz));
            WL_TZ(kk,ii)=WLz(Tz{jj}(kk,Bz(ii):Dz));
            
        end
        ZC{jj}=ZC_TZ;
        MAV{jj}=MAV_TZ;
        SSc{jj}=SSC_TZ;
        WL{jj}=WL_TZ;
        
    end
    
end
cord=0;
%% covariance
for ii=1:numberOfFiles
    clear a
    clear V
    clear D
    a=[WL{ii};SSc{ii};MAV{ii};ZC{ii}];
    cord=cord+1;
    cov_mat{cord}=cov(a');
    a_mat{ii}=a;
    [V,D] = eig(cov_mat{ii});
    eig_vec{ii}=V;
    eig_val{ii}=D;
end
eig_val2=eig_val;%just so we don't mess up anything
for ii=1:numberOfFiles
    kk=0;
    eig_cut=10^-6;
    for jj=1:32
        if eig_val2{ii}(jj,jj)>eig_cut
            eig_highest=eig_val2{ii}(jj,jj);
            remove=jj;
            kk=kk+1;
            eig_val2{ii}(remove,remove)=0;
            eig_new{ii}(kk,:)=eig_val{ii}(remove,remove);
        end
        
    end
    
    
end

%% can someone start the Euclidean? This one is stumping me

