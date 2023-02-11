function [XN,YN]=SmothBound_N2(XS,YS,N)

    DS=sqrt((XS(2:end)-XS(1:(end-1))).^2+(YS(2:end)-YS(1:(end-1))).^2);
    CS=[0 cumsum(DS)];
    L=sum(DS); 
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

    %XN=fliplr(XN); YN=fliplr(YN);
    