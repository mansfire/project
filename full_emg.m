clc;
clear all;
    close all;
%% hand close%%
%% Training data%%
load('Hand Close Training.mat')
fs=samplingRate ;%sampling rate
dt=1/fs;
time=(dt:dt:length_sec);
mymodel.name{1}='hand_close_train';
mymodel.data{1}=Data;
%% Validation
load('Hand Close Validation.mat')
mymodel.name{2}='hand_close_val';
mymodel.data{2}=Data;
%% hand open
%% Training data%%
load('Hand Open Train.mat')
mymodel.name{3}='hand_open_train';
mymodel.data{3}=Data;
%% Validation
load('Hand Open Validation.mat')
mymodel.name{4}='hand_open_val';
mymodel.data{4}=Data;
%% off
%% Training data%%
load('Off Train.mat')
mymodel.name{5}='off_train ';
mymodel.data{5}=Data;
%% Validation
load('Off Validation.mat')
mymodel.name{6}='off_val';
mymodel.data{6}=Data;
%% thumb ab
%% Training data%%
load('Thumb Abduction Training.mat')
mymodel.off_train=Data;
mymodel.name{7}='thumb_ab_train';
mymodel.data{7}=Data;
%% Validation
load('Thumb Abduction Validation.mat')
mymodel.name{8}='thumb_ab_val';
mymodel.data{8}=Data;
%% thumb add
%% Training data%%
load('Thumb Adduction Training.mat')
mymodel.name{9}='thumb_add_train';
mymodel.data{9}=Data;
%% Validation
load('Thumb Adduction Validation.mat')
mymodel.name{10}='thumb_add_val';
mymodel.data{10}=Data;
%% wrist ext
%% Training data%%
load('Wrist Extension Training.mat')
mymodel.name{11}='wrist_ext_train';
mymodel.data{11}=Data;
%% Validation
load('Wrist Extension Validation.mat')
mymodel.name{12}='wrist_ext_val';
mymodel.data{12}=Data;
%% wrist flex
%% Training data%%
load('Wrist Flexion Training.mat')
mymodel.name{13}='wrist_flex_train';
mymodel.data{13}=Data;
%% Validation
load('Wrist Flexion Validation.mat')
mymodel.name{14}='wrist_flex_val';
mymodel.data{14}=Data;
%% wrist pronation
%% Training data%%
load('Wrist Pronation Training.mat')
mymodel.name{15}='wrist_pronation_train';
mymodel.data{15}=Data;
%% Validation
load('Wrist Pronation Validation.mat')
mymodel.name{16}='wrist_pronation_val';
mymodel.data{16}=Data;
%% wrist radial dev
%% Training data%%
load('Wrist Rad Dev Training.mat')
mymodel.name{17}='wrist_radial_dev_train';
mymodel.data{17}=Data;
%% Validation
load('Wrist Rad Dev Validation.mat')
mymodel.name{18}='wrist_radial_dev_val';
mymodel.data{18}=Data;
%% wrist supination
%% Training data%%
load('Wrist Supination Training.mat')
mymodel.name{19}='wrist_supination_train';
mymodel.data{19}=Data;
%% Validation
load('Wrist Supination Validation.mat')
mymodel.name{20}='wrist_supination_val';
mymodel.data{20}=Data;
%% wrist ulnar dev
%% Training data%%
load('Wrist Ulnar Dev Training.mat')
mymodel.name{21}='wrist_ulnar_dev_train';
mymodel.data{21}=Data;
%% Validation
load('Wrist Ulnar Dev Validation.mat')
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
     data=[data{1}, data{2}, data{3},data{4},data{5},data{6},data{7},data{8}];
clear filtdata1;
clear filtdata2;
clear filtdata3;
clear filtdata4;
clear fulldata;
bp1  = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp, ...
         'SampleRate',fs);
%% band pass

% Filter the data and compensate for delay
D1 = round(mean(grpdelay(bp1))); % filter delay
a=0;

    bpfilt = filter(bp1,[data; zeros(D1,8)]);
    filtdata1 = bpfilt(D1+1:end,:);
    

n1f1=n1-1;
n1f2=n1+1;
notch1 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n1f1,'PassbandFrequency2',n1f2,'SampleRate',fs);

%% notch 60hz
% Filter the data and compensate for delay
D2 = round(mean(grpdelay(notch1))); % filter delay


    bpnotch = filter(notch1,[filtdata1(:,:); zeros(D2,8)]);
    filtdata2 = bpnotch(D2+1:end,:);


n2f1=n2-1;
n2f2=n2+1;
notch2 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n2f1,'PassbandFrequency2',n2f2,'SampleRate',fs);

%% notch 120hz
% Filter the data and compensate for delay
D3 = round(mean(grpdelay(notch2))); % filter delay

    bpnotch2 = filter(notch2,[filtdata2(:,:); zeros(D3,8)]);
    filtdata3 = bpnotch2(D3+1:end,:);


n3f1=n3-1;
n3f2=n3+1;
notch3 = designfilt('bandstopiir', ...
  'FilterOrder',16, ...
  'PassbandFrequency1',n3f1,'PassbandFrequency2',n3f2,'SampleRate',fs);

%% notch 180hz
% Filter the data and compensate for delay
D4 = round(mean(grpdelay(notch3))); % filter delay

    bpnotch3 = filter(notch3,[filtdata3(:,:); zeros(D4,8)]);
    filtdata4 = bpnotch3(D4+1:end,:);

fulldata=filtdata4;

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
total_data=[];
for ii=1:22
    total=zeros(1,length(ylp.data{ii}(:,1)));
    for jj=1:8
        for kk=1:length(ylp.data{ii}(:,jj))
            total(:,kk)=total(:,kk)+abs(ylp.data{ii}(kk,jj));
        end
    end
    total_data{ii}=total;
end


for ii=1:22
    Smoothed_data{ii}=smoothdata(total_data{ii},'gaussian',150);
    std_data(ii)=std(Smoothed_data{ii});
    mean_data(ii)=mean(Smoothed_data{ii});
end
S=0;
j=0;
Start=[];
Start1=[];
Stop=[];
Stop1=[];
data1=smoothdata(abs(ylp.data{1}(:,2)),'gaussian',300);
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
data2=smoothdata(abs(ylp.data{2}(:,2)),'gaussian',300);
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
data3=smoothdata(abs(ylp.data{3}(:,5)),'gaussian',300);
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
data4=smoothdata(abs(ylp.data{4}(:,5)),'gaussian',300);
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
data7=smoothdata(abs(ylp.data{7}(:,8)),'gaussian',900);
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
data8=smoothdata(abs(ylp.data{8}(:,8)),'gaussian',3000);
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
filename='smoothed.mat';
save(filename)
%% new data set


for ii=1:22
    k=1;
    start=Start{ii};
    stop=Stop{ii};
    for jj=1:length(start)
        
        for kk=start(jj):stop(jj)
            New{ii}(k)=Smoothed_data{ii}(kk);
            k=k+1;
        end
    end
end


        
%% index on data
p=[];
for jj=1:22
    clear p;
    start=Start{ii};
    stop=Stop{ii};
    for ii=1:length(Start1)
        p=[p Start1(ii):Stop1(ii)];
    end
    P{ii}=p;
end
%% create trimmed data set
TrimmedTF=[];
 for ii=1:22
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
 T=zeros(8,240100);
 for kk=1:22
     for ii=1:8
         for jj=1:length(ylp.data{kk}{ii})
             T(ii,jj)=(ylp.data{kk}{ii}(jj));
         end
     end
     Tz{kk}=T;
 end
 for jj=1:22
     for kk=1:8
     nBin=floor(length(Tz{jj})/binsize);
     Bz=floor(linspace(1,length(Tz{jj}),nBin));
     Ez=[];
     ZC_TZ=Ez;
     SSC_TZ=Ez;
     MAV_TZ=Ez;
     WL_TZ=Ez;
     for ii=1:nBin-1
         Dz=Bz(ii+1);
        
         ZC_TZ(kk,ii)=ZCz(Tz{jj}(kk,Bz(ii):Dz));
         MAV_TZ(kk,ii)=MAVz(Tz{jj}(kk,Bz(ii):Dz));
         SSC_TZ(kk,ii)=SSCz(Tz{jj}(kk,Bz(ii):Dz));
         WL_TZ(kk,ii)=WLz(Tz{jj}(kk,Bz(ii):Dz));

     end
     ZC{jj}=ZC_TZ(:,1:1596);
     MAV{jj}=MAV_TZ(:,1:1596);
     SSc{jj}=SSC_TZ(:,1:1596);
     WL{jj}=WL_TZ(:,1:1596);
     end
 end
 cord=0;
%% covariance
for ii=1:22
    clear a
    clear V
    clear D
    a=[WL{ii};SSc{ii};MAV{ii};ZC{ii}];
    cord=cord+1;
    cov_mat{cord}=cov(a');   
    [V,D] = eig(cov_mat{ii});
    eig_vec{ii}=V;
    eig_val{ii}=D;
end
for ii=1:22
    mean_class=mean(V*a);
    dif=0;
    for jj=1:32
        for kk=1:length(a(ii))
            dif=dif+(abs(mean_class(kk)-a(jj,kk)))^2;
        end
    end
end
