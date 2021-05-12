% Feature extraction function
% Calculates mean absolute value of given vector of data

function mav_out=MAVz(x)
L=length(x);
tot_out=0;
for ii=1:length(x)
    tot_out=tot_out+abs(x(ii));
end
mav_out=tot_out/L;
end