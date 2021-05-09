S=0;
j=0;
Start=[];
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
