% Feature extraction function
% Calculates waveform length of given vector of data

function wl_out=WLz(x)
wl_out=0;
for ii=1:length(x)-1
    dxk=abs(x(ii+1)-x(ii));
    wl_out=wl_out+dxk;
end
end