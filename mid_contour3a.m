function [Xm,Ym,Lm] = mid_contour3a(X1,Y1,X2,Y2,str)

    Xm = -ones(size(X1));
    Ym = Xm; 
    L1 = length(X1);
    L2 = length(X2);

    wb = waitbar(0,str);
    for i = 1:L1
        waitbar(i/L1);
        
        dis = (X1(i)-X2).^2 + (Y1(i)-Y2).^2;
        
        [~,si]  = sort(dis);
        X2 = X2(si);
        Y2 = Y2(si);
        
        
        %crs = zeros(1,L2);
        %done = 0;
        for m = 1:L2    
            
            q1 = [X1(i), Y1(i)]; 
            q2 = [X2(m), Y2(m)];
            a2x = q1(1) - q2(1);
            a2y = q1(2) - q2(2);
            
            x1 = X1; 
            y1 = Y1;
            x2 = [X1(2:end), X1(1)];
            y2 = [Y1(2:end), Y1(1)];
            
            a1x = x1 - q2(1);
            a1y = y1 - q2(2);
            s1 = a1x*a2y - a1y*a2x;

            a1x = x2 - q2(1);
            a1y = y2 - q2(2);
            s2 = a1x*a2y - a1y*a2x;
            
            s = s1.*s2;
            
            IND = find(s<=0);
            if i == 1
                IND = IND(IND~=1 & IND~=L1);
            else
                IND = IND(IND~=i & IND~=(i-1));
            end
            
            chk = zeros(size(IND));
            %xc = zeros(size(IND));
            %yc = zeros(size(IND));
            for k = 1:length(IND)
                r1 = [x1(IND(k)) y1(IND(k))]; r2 = [x2(IND(k)) y2(IND(k))];
                [~,~,chk(k)] = intercheck(q1,q2,r1,r2);
                %[xc(k),yc(k),chk(k)] = intercheck(q1,q2,r1,r2);
                
            end
            if sum(chk)==0
                Xm(i) = (X1(i)+X2(m))/2;
                Ym(i) = (Y1(i)+Y2(m))/2;
                %done = 1;
                %plot([X1(i) X2(m)],[Y1(i) Y2(m)],'Color',[.5 .5 .5]);
                
                break;
            end
           
        end
        
        %if i>1 && done == 0
        %    Xm(i) = Xm(i-1);
        %    Ym(i) = Ym(i-1);
        %elseif i == 1 && done == 0 
        %    needtofix
        %end
            
        
        %dis_non = dis(crs==0);
        %X2_non = X2(crs==0);
        %Y2_non = Y2(crs==0);
        
        %[~,jj] = min(dis_non);
        %if isempty(jj) && i>1
        %    Xm(i) = Xm(i-1);
        %    Ym(i) = Ym(i-1);
        %else
        %    Xm(i) = (X1(i)+X2_non(jj))/2;
        %    Ym(i) = (Y1(i)+Y2_non(jj))/2;
        %end
        
    end

    Xm = Xm(Xm>=0);
    Ym = Ym(Ym>=0);
    
    
    x = [Xm Xm(1)];
    y = [Ym Ym(1)];
    
    d=sqrt((x(2:end)-x(1:(end-1))).^2+(y(2:end)-y(1:(end-1))).^2);
    
    x = x(d>1e-6);
    y = y(d>1e-6);
    
    
    Lm = sum(d);
    n = ceil(Lm);
    [Xm,Ym]=SmothBound_N2(x,y,n);
    
    close(wb);

    