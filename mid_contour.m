function [Xm,Ym,Lm] = mid_contour(X1,Y1,X2,Y2)

    Xm = zeros(size(X1));
    Ym = Xm;

    for i = 1:length(X1)
        dis = (X1(i)-X2).^2 + (Y1(i)-Y2).^2;
        [~,j] = min(dis);
        Xm(i) = (X1(i)+X2(j))/2;
        Ym(i) = (Y1(i)+Y2(j))/2;
    end

    x = [Xm Xm(1)];
    y = [Ym Ym(1)];
    d=sqrt((x(2:end)-x(1:(end-1))).^2+(y(2:end)-y(1:(end-1))).^2);
    Lm = sum(d);
    n = ceil(Lm);
    [Xm,Ym]=SmothBound_N2(x,y,n);

