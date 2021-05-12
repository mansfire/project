% Feature extraction function
% Calculates number of slope sign change of a given vector of data

function ssc_out=SSCz(x)
e=1*10^-5;
ssc_out=0;
for ii=2:length(x)-1
    slope1=abs(x(ii)-x(ii-1));
    slope2=abs(x(ii+1)-x(ii));
    if ((x(ii)>x(ii+1)&&x(ii)>x(ii-1))||(x(ii)<x(ii+1)&&x(ii)<x(ii-1))) && (slope1>=e||slope2>=e)
        ssc_out=ssc_out+1;
    end
end
end