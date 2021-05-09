function ZC_out=ZCz(x)
e=1*10^-5;
if x(1)==0
    ZC_out=1;
else
    ZC_out=0;
end
z_prev=abs(x(1));
for ii=2:length(x)
    z_cur=abs(x(ii));
    z_dif=abs(z_cur-z_prev);
    if z_cur<=e &&z_prev>e && z_dif>=e
        ZC_out=ZC_out+1;
    end
end
end
