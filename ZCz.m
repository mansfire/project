% Feature extraction function
% Calculates number of zero crossings for given vector of data

function ZC_out=ZCz(x)
e=1*10^-5;
if x(1)==0
    ZC_out=1;
else
    ZC_out=0;
end
z_prev=(x(1));
for ii=2:length(x)
    z_cur=(x(ii));
    z_dif=abs(z_cur-z_prev);
    if ((z_cur<=0 &&z_prev>0)|| (z_cur>0 &&z_prev<=0)) && z_dif>=e
        ZC_out=ZC_out+1;
    end
end
end
