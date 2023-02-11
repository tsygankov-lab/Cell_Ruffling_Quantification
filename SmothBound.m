function [XN,YN,XS,YS]=SmothBound(F,K)

    f=figure;

    F=im2bw(F,0.5); F=+F;
    h=fspecial('gaussian', [5 5], 1);
    G=imfilter(F, h,'replicate');
    SmBound = contour(G,1);  
    
    ii=find(SmBound(1,:)<1);
    [vv, jj]=max(SmBound(2,ii));
    
    bg=ii(jj)+1;
    en=ii(jj)+vv;
    
    XS=SmBound(1,bg:en); YS=SmBound(2,bg:en); %NS=SmBound(2,1);

    DS=sqrt((XS(2:end)-XS(1:(end-1))).^2+(YS(2:end)-YS(1:(end-1))).^2);
    CS=[0 cumsum(DS)];
    L=sum(DS); N=K*round(L);
    r=L/N;

    XN=zeros(1,N); YN=zeros(1,N);
    RD=0; XN(1)=XS(1); YN(1)=YS(1);
    for i=1:(N-1)    
        RD=RD+r;
        [vi, mi]=max(CS(CS<=RD));
        w=RD-vi;
        XN(i+1)=XS(mi)+(XS(mi+1)-XS(mi))*w/DS(mi);
        YN(i+1)=YS(mi)+(YS(mi+1)-YS(mi))*w/DS(mi);    
    end

    XN=fliplr(XN); YN=fliplr(YN);
    close(f);