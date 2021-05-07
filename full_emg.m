
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
    %mymodel.data{idx}{6} = zeros(length(mymodel.data{idx}{6}), 1); %
    %channel 6 might be bad. Remove?
    mymodel.length{idx} = emgDataStruct.length_sec;
end

%% Construct the filters
lp=480;
hp=30;

bp1  = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp, ...
    'SampleRate',fs);

total_data = cell(1,numberOfFiles);
for kk=1:numberOfFiles
    [~, name, ~] = fileparts(fdsGroup4.Files{kk});
    figure;
    data=mymodel.data{kk};
    data=cell2mat(data); % converts data into matrix format.
    
    %% band pass
    fulldata = filtfilt(bp1, data); % No phase shift.
    
    % Spectrum interpolation to remove 60 hz + harmonic noise. Works in the
    % frequency domain and does not rely on filters. No phase shift.
    for idx = 1:width(data)
        fulldata(:, idx) = spectrumInterpolation(fulldata(:, idx), fs, 60, 3, 2);
        nexttile;
        plot(fulldata(:,idx));
        title(name);
    end
    mymodel.data{kk}=fulldata; % write back to conserve memory.
    mymodel.meanChan{kk} = mean(fulldata,2);
    nexttile;
    plot(mymodel.meanChan{kk});
    title(name);
end

%% add all the channels for each of the 11 test and validation postures
for ii=1:numberOfFiles
    total_data{idx}=sum(abs(mymodel.data{idx}), 2); % Rectifies the data, then adds all channels together.
end

%% smooth out the data to make it easuer to analyze
Smoothed_data = cell(1,numberOfFiles);
for ii=1:numberOfFiles
    Smoothed_data{ii}=smoothdata(total_data{ii},'gaussian',150);
    std_data(ii)=std(Smoothed_data{ii});
    mean_data(ii)=mean(Smoothed_data{ii});
end
%% for each posture, try to capture apptox. 50% (0.5) of the data, us whichever channel is "prettiest"
S=0;
Start1=[];
Stop=[];
Stop1=[];
j=0;
Start=[];
data1=smoothdata(abs(mymodel.data{1}(:,2)),'gaussian',300);
mean1=mean(data1);
std1=std(data1);
for jj=300:length(data1)
    if (data1(jj)>1.5*mean1+1.5*std1) &&(S==0)
        S=1;
        j=j+1;
        Start1=[Start1 jj];
    elseif (data1(jj)<mean1-1.1*std1) &&(S==1)
        S=0;
        Stop1=[Stop1 jj];
        
    end
    
end
if length(Stop1)<length(Start1)
    Stop1=[Stop1 240100];
end
Start{1}=Start1;
Stop{1}=Stop1;
S=0;
j=0;

Start2=[];

Stop2=[];
data2=smoothdata(abs(mymodel.data{2}(:,2)),'gaussian',300);
mean2=mean(data2);
std2=std(data2);
for jj=300:length(data2)
    if (data2(jj)>1.5*mean2+std2) &&(S==0)
        S=1;
        j=j+1;
        Start2=[Start2 jj];
    elseif (data2(jj)<mean2-1.1*std2) &&(S==1)
        S=0;
        Stop2=[Stop2 jj];
        
    end
    
end
if length(Stop2)<length(Start2)
    Stop2=[Stop2 jj];
end

Start{2}=Start2;
Stop{2}=Stop2;

S=0;
j=0;

Start3=[];

Stop3=[];
data3=smoothdata(abs(mymodel.data{3}(:,5)),'gaussian',300);
mean3=mean(data3);
std3=std(data3);
for jj=300:length(data3)
    if (data3(jj)>1.7*mean3+1.3*std3) &&(S==0)
        S=1;
        j=j+1;
        Start3=[Start3 jj];
    elseif (data3(jj)<mean3-0.9*std3) &&(S==1)
        S=0;
        Stop3=[Stop3 jj];
        
    end
    
end
if length(Stop3)<length(Start3)
    Stop3=[Stop3 jj];
end
Start{3}=Start3;
Stop{3}=Stop3;

S=0;
j=0;

Start4=[];

Stop4=[];
data4=smoothdata(abs(mymodel.data{4}(:,5)),'gaussian',300);
mean4=mean(data4);
std4=std(data4);
for jj=300:length(data4)
    if (data4(jj)>1.3*mean4+1.6*std4) &&(S==0)
        S=1;
        j=j+1;
        Start4=[Start4 jj];
    elseif (data4(jj)<mean4-std4) &&(S==1)
        S=0;
        Stop4=[Stop4 jj];
        
    end
    
end
Start{4}=Start4;
Stop{4}=Stop4;

S=0;
j=0;

Start5=[300];

Stop5=[249399];


Start{5}=Start5;
Stop{5}=Stop5;
Start6=[300];

Stop6=[249399];


Start{6}=Start6;
Stop{6}=Stop6;
S=0;
j=0;

Start7=[];

Stop7=[];
data7=smoothdata(abs(mymodel.data{7}(:,8)),'gaussian',900);
mean7=mean(data7);
std7=std(data7);
for jj=300:length(data7)
    if (data7(jj)>1.1*mean7+1.2*std7) &&(S==0)
        S=1;
        j=j+1;
        Start7=[Start7 jj];
    elseif (data7(jj)<(mean7-std7)) &&(S==1)
        S=0;
        Stop7=[Stop7 jj];
        
    end
    
end
if length(Stop7)<length(Start7)
    Stop7=[Stop7 jj];
end
Start{7}=Start7;
Stop{7}=Stop7;


S=0;
j=0;

Start8=[];

Stop8=[];
data8=smoothdata(abs(mymodel.data{8}(:,8)),'gaussian',3000);
mean8=mean(data8);
std8=std(data8);
for jj=300:length(data8)
    if (data8(jj)>mean8+0.6*std8) &&(S==0)
        S=1;
        j=j+1;
        Start8=[Start8 jj];
    elseif (data8(jj)<mean8) &&(S==1)
        S=0;
        Stop8=[Stop8 jj];
        
    end
    
end
if length(Stop8)<length(Start8)
    Stop8=[Stop8 jj];
end
Start{8}=Start8;
Stop{8}=Stop8;

S=0;
j=0;

Start9=[];

Stop9=[];
data9=smoothdata(abs(mymodel.data{9}(:,4)),'gaussian',5000);
mean9=mean(data9);
std9=std(data9);
for jj=300:length(data9)
    if (data9(jj)>mean9) &&(S==0)
        S=1;
        j=j+1;
        Start9=[Start9 jj];
    elseif (data9(jj)<mean9-0.33*std9) &&(S==1)
        S=0;
        Stop9=[Stop9 jj];
        
    end
    
end
if length(Stop9)<length(Start9)
    Stop9=[Stop9 jj];
end
Start{9}=Start9;
Stop{9}=Stop9;
S=0;
j=0;

Start10=[];

Stop10=[];
data10=smoothdata(abs(mymodel.data{10}(:,4)),'gaussian',5000);
mean10=mean(data10);
std10=std(data10);
for jj=300:length(data10)
    if (data10(jj)>mean10) &&(S==0)
        S=1;
        j=j+1;
        Start10=[Start10 jj];
    elseif (data10(jj)<mean10-0.43*std10) &&(S==1)
        S=0;
        Stop10=[Stop10 jj];
        
    end
    
end
if length(Stop10)<length(Start10)
    Stop10=[Stop10 jj];
end
Start{10}=Start10;
Stop{10}=Stop10;

filename='smoothed.mat';
save(filename)
%% new data set


for ii=1:numberOfFiles
    k=1;
    start=Start{ii};
    stop=Stop{ii};
    
    for jj=1:length(start)%look at each start point
        
        for kk=start(jj):stop(jj)%look at the data between start 1 and stop 1, etc.
            
            New{ii}(k)=Smoothed_data{ii}(kk);
            k=k+1;
        end
    end
end



%% index on data
p=[];
for jj=1:numberOfFiles
    clear p;
    start=Start{ii};
    stop=Stop{ii};
    for ii=1:length(Start1)
        p=[p Start1(ii):Stop1(ii)];
    end
    P{ii}=p;
end
%% create trimmed data set% this looks at numberOfFiles postures (11 X2) and 8 channels and only takes th on aspects of each
TrimmedTF=[];
for ii=1:numberOfFiles
    for jj=1:8
        for kk=1:length(p)
            Trim(jj,kk)=mymodel.data{ii}{jj}(P{ii}(kk));
        end
    end
    TrimmedTF{ii}=Trim;
end

filename='trimmed.mat';
save(filename)
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

