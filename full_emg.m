clc;
clear all;
close all;
%% hand close%%
%% Training data%%
load('hand_close_train.mat')
fs=samplingRate ;%sampling rate
dt=1/fs;
time=(dt:dt:length_sec);
mymodel.name{1}='hand_close_train';
mymodel.data{1}=Data;
%% Validation
load('hand_close_val.mat')
mymodel.name{2}='hand_close_val';
mymodel.data{2}=Data;
%% hand open
%% Training data%%
load('hand_open_train.mat')
mymodel.name{3}='hand_open_train';
mymodel.data{3}=Data;
%% Validation
load('hand_open_val.mat')
mymodel.name{4}='hand_open_val';
mymodel.data{4}=Data;
%% off
%% Training data%%
load('off_train.mat')
mymodel.name{5}='off_train ';
mymodel.data{5}=Data;
%% Validation
load('off_val.mat')
mymodel.name{6}='off_val';
mymodel.data{6}=Data;
%% thumb ab
%% Training data%%
load('thumb_ab_train.mat')
mymodel.off_train=Data;
mymodel.name{7}='thumb_ab_train';
mymodel.data{7}=Data;
%% Validation
load('thumb_ab_val.mat')
mymodel.name{8}='thumb_ab_val';
mymodel.data{8}=Data;
%% thumb add
%% Training data%%
load('thumb_add_train.mat')
mymodel.name{9}='thumb_add';
mymodel.data{9}=Data;
%% Validation
load('thumb_add_val.mat')
mymodel.name{10}='thumb_add_val';
mymodel.data{10}=Data;
%% wrist ext
%% Training data%%
load('wrist_ext_train.mat')
mymodel.name{11}='wrist_ext_train';
mymodel.data{11}=Data;
%% Validation
load('wrist_ext_val.mat')
mymodel.name{12}='wrist_ext_val';
mymodel.data{12}=Data;
%% wrist flex
%% Training data%%
load('wrist_flex_train.mat')
mymodel.name{13}='wrist_flex_train';
mymodel.data{13}=Data;
%% Validation
load('wrist_flex_val.mat')
mymodel.name{14}='wrist_flex_val';
mymodel.data{14}=Data;
%% wrist pronation
%% Training data%%
load('wrist_pronation_train.mat')
mymodel.name{15}='wrist_pronation_train';
mymodel.data{15}=Data;
%% Validation
load('wrist_pronation_val.mat')
mymodel.name{16}='wrist_pronation_val';
mymodel.data{16}=Data;
%% wrist radial dev
%% Training data%%
load('wrist_radial_dev_train.mat')
mymodel.name{17}='wrist_radial_dev_train';
mymodel.data{17}=Data;
%% Validation
load('wrist_radial_dev_val.mat')
mymodel.name{18}='wrist_radial_dev_val';
mymodel.data{18}=Data;
%% wrist supination
%% Training data%%
load('wrist_supination_train.mat')
mymodel.name{19}='wrist_supination_train';
mymodel.data{19}=Data;
%% Validation
load('wrist_supination_val.mat')
mymodel.name{20}='wrist_supination_val';
mymodel.data{20}=Data;
%% wrist ulnar dev
%% Training data%%
load('wrist_ulnar_dev_train.mat')
mymodel.name{21}='wrist_ulnar_dev_train';
mymodel.data{21}=Data;
%% Validation
load('wrist_ulnar_dev_val.mat')
mymodel.name{22}='wrist_ulnar_dev_val.mat';
mymodel.data{22}=Data;

%% filters
lp=480;
hp=30;
n1=60;
n2=120;
n3=180;
for kk=1:22
    data=mymodel.data{kk};
for ll=1:8
bp1  = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp, ...
         'SampleRate',fs);
%% band pass
[P,F] = pwelch(data{ll},ones(4097,1),round(4097/2),4097,fs,'power');

% Filter the data and compensate for delay
D1 = round(mean(grpdelay(bp1))); % filter delay
a=0;
for jj=1:8
    bpfilt = filter(bp1,[Data{jj}; zeros(D1,1)]);
    filtdata1(:,jj) = bpfilt(D1+1:end);
    
end
n1f1=n1-1;
n1f2=n1+1;
notch1 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n1f1,'PassbandFrequency2',n1f2,'SampleRate',fs);

%% notch 60hz
% Filter the data and compensate for delay
D2 = round(mean(grpdelay(notch1))); % filter delay

for jj=1:8
    bpnotch = filter(notch1,[filtdata1(:,jj); zeros(D2,1)]);
    filtdata2(:,jj) = bpnotch(D2+1:end);
end

n2f1=n2-1;
n2f2=n2+1;
notch2 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n2f1,'PassbandFrequency2',n2f2,'SampleRate',fs);

%% notch 120hz
% Filter the data and compensate for delay
D3 = round(mean(grpdelay(notch2))); % filter delay

for jj=1:8
    bpnotch2 = filter(notch2,[filtdata2(:,jj); zeros(D3,1)]);
    filtdata3(:,jj) = bpnotch2(D3+1:end);
end

n3f1=n3-1;
n3f2=n3+1;
notch3 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n3f1,'PassbandFrequency2',n3f2,'SampleRate',fs);

%% notch 180hz
% Filter the data and compensate for delay
D4 = round(mean(grpdelay(notch3))); % filter delay

for jj=1:8
    bpnotch3 = filter(notch3,[filtdata3(:,jj); zeros(D4,1)]);
    filtdata4(:,jj) = bpnotch3(D4+1:end);
end
fulldata{:,ll}=filtdata4;
end
ylp.data{kk}=fulldata;
end
filename='filtered.mat';
save(filename)
% for ii=1:12
%     data=ylp.data{ii};
%     for jj=1:8
%         mydata=(abs(data{jj})./4) ;
%         bindata{:,jj}=mydata;
%     end
%     ylp2.data{jj}=bindata;
% end
% for jj=1:12
% binned=ylp.data{jj};
% data_bin=binned{1};
% figure(jj)
% for ii = 1: 4
% subplot(4,2,ii);
% plot(data_bin(:,ii))
%     if ii == 1
%         ylabel('Voltage [mV]')
%     end
%     ylim([0, 200])
%     title(['Bin EMG ', num2str(ii), 'for ', mymodel.name{jj}])
% end
% hold on
% for ii = 1: 4
%     subplot(4,2,ii+4);
%     plot(data_bin(:,ii+4))
%     ylim([0, 200])
%     title(['Bin EMG ', num2str(ii+4),'for ' mymodel.name{jj}])
% end
% hold off
% end
total_data=zeros(22,240100);
for ii=1:22
    for jj=1:8
        for kk=1:240100
            total_data(ii,kk)=total_data(ii,kk)+abs(mymodel.data{ii}{jj}(kk));
        end
    end
end


for ii=1:22
    Smoothed_data(ii,:)=smoothdata(total_data(ii,:),'gaussian',150);
    std_data(ii)=std(Smoothed_data(ii,1:240100));
    mean_data(ii)=mean(Smoothed_data(ii,1:240100));
end
S=0;
j=0;
Start=[];
Start1=[];
Stop=[];
Stop1=[];
for jj=300:240100
   if (Smoothed_data(1,jj)>1.03*mean_data(1)+std_data(1))&&(S==0)
         S=1;
         j=j+1;
         Start1=[Start1 jj];
     elseif (Smoothed_data(1,jj)<0.97*mean_data(1)- std_data(1))&&(S==1)
         S=0;
         Stop1=[Stop1 jj];
            
   end

end
Start{1}=Start1;
Stop{1}=Stop1;
Start2=[];
Stop2=[];
S=0;
j=0;
for jj=300:240100
   if (Smoothed_data(2,jj)>1.03*mean_data(2)+1.9*std_data(2))&&(S==0)
         S=1;
         j=j+1;
         Start2=[Start2 jj];
     elseif (Smoothed_data(2,jj)<0.96*mean_data(2)- std_data(2))&&(S==1)
         S=0;
         Stop2=[Stop2 jj];
            
    end
end
Start{2}=Start2;
Stop{2}=Stop2;
Start3=[];
Stop3=[];
S=0;
j=0;
for jj=300:240100
   if (Smoothed_data(3,jj)>1.07*mean_data(3)+std_data(2))&&(S==0)
         S=1;
         j=j+1;
         Start3=[Start3 jj];
     elseif (Smoothed_data(3,jj)<mean_data(3)- 2*std_data(3))&&(S==1)
         S=0;
         Stop3=[Stop3 jj];
            
    end
end
Start{3}=Start3;
Stop{3}=Stop3;
filename='smoothed.mat';
save(filename)
%% new data set

New=zeros(22,819);
for ii=1:22
    
    for jj=1:length(Start2)
        k=1;
        for kk=Start2(jj):Stop2(jj)
            New(ii,k)=Smoothed_data(ii,kk);
            k=k+1;
        end
    end
end


        
%% index on data
p=[];
for ii=1:length(Start2)
    p=[p Start2(ii):Stop2(ii)];
end
%% create trimmed data set
TrimmedTF=[];
 for ii=1:22
     for jj=1:8
         for kk=1:length(p)
             Trim(jj,kk)=mymodel.data{ii}{jj}(p(kk));
         end
     end
     TrimmedTF(ii,:,:)=Trim;
 end
 
 filename='trimmed.mat';
save(filename)
 %% off data prep
 binsize=0.05*fs;
 T=zeros(8,length(mymodel.data{9}{1}));
 for kk=1:22
     for ii=1:8
         for jj=1:240100
             T(ii,jj)=(mymodel.data{kk}{ii}(jj));
         end
     end
     Tz(kk,:,:)=T;
 end
 for jj=1:22
     nBin=floor(length(Tz(jj,:,:))/binsize);
     Bz=floor(linspace(1,length(Tz(jj,:,:)),nBin));
     Ez=[];
     ZC_TZ=Ez;
     SSC_TZ=Ez;
     MAV_TZ=Ez;
     WL_TZ=Ez;
     for ii=1:nBin-1
         Dz=Bz(ii+1);
        
         ZC_TZ(:,ii)=ZCz(Tz(jj,:,Dz));
         MAV_TZ(:,ii)=MAVz(Tz(jj,:,Dz));
         SSC_TZ(:,ii)=SSCz(Tz(jj,:,Dz));
         WL_TZ(:,ii)=WLz(Tz(jj,:,Dz));

     end
     ZC(jj,:,:)=ZC_TZ;
     MAV(jj,:,:)=MAV_TZ;
     SSc(jj,:,:)=SSC_TZ;
     WL(jj,:,:)=WL_TZ;
 end
 
 
