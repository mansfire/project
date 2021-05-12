% LDA
% TODO: covariance only on training data (todo)
cord=0;
for ii=1:numberOfPoses
    a=[WL_train{ii};SSC_train{ii};MAV_train{ii};ZC_train{ii}];%create a matrix of our feature
    cord=cord+1; %%increment the index
    [M,N]=size(a);
    mn=mean(a,2);
    a=a-repmat(mn,1,N);
    cov_mat{cord}=cov(a*a');%%create the covariance natrix
    a_mat{ii}=a;%%save feature

end

for ii=1:numberOfPoses
% TODO: W matrix only on training data (todo)
%     a_mat{ii} = abs(a_mat{ii});
    u=mean(a_mat{ii}');
    class_means{ii} = mean(a_mat{ii},2);
    bcs=cov(u'*u);%%between class
    wcs=zeros(32,32);%%lets do within class now
    for jj=1:numberOfPoses
        wcs=wcs+cov_mat{jj};
    end
    wcs=wcs/(numberOfPoses-1);
    Op=wcs'*bcs;%%optimization matrix
    [Z,W]=eig(Op);
    W=real(W);
    Z=real(Z);
 
    Wact=diag(W);%%look at the acual values
    [val,indexT]=sort(-1*Wact);

    Windex=Wact(indexT);
    Z=Z(:,indexT);
    eig_vec{ii}=Z;%%eigen vector

    Y{ii}=Z'*a_mat{ii};%%use this to find euclidean
    Y_avg{ii} = mean(Y{ii},2);
end