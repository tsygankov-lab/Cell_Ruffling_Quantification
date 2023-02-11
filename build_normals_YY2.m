function [elapsedTime,N_noINT,convg] = build_normals_YY2(fn1,fn2,dp,maxitr,save_fold)

    tic
    
    file_location1 = fn1;
    file_location2 = fn2;

    info = imfinfo(file_location1);
    T = 20;
    cont_N = fix(length(info)/T);
    if mod(length(info),T) == 0
        cont_N = cont_N - 1;
    end
    se = strel('disk',1);

    %dp = 1;
    convg = zeros(1,cont_N);
    N_noINT = zeros(1,cont_N);


    cellXin = cell(1,cont_N);
    cellYin = cell(1,cont_N);
    cellXou = cell(1,cont_N);
    cellYou = cell(1,cont_N);
    
    cellQ1 = cell(1,cont_N);
    cellQ2 = cell(1,cont_N);
    cellXs = cell(1,cont_N);
    cellYs = cell(1,cont_N);
    cellDs = cell(1,cont_N);
    cellXmd = cell(1,cont_N);
    cellYmd = cell(1,cont_N);
    noINTER = cell(1,cont_N);

    for uu = 1:cont_N
        %% Finding outer, inner, first, and last boundaries

        st = 1+T*(uu-1);
        en = 1+T*uu;

        wb = waitbar(0,['time block ' num2str(uu) ', boundaries']);
        for fr = st:en

            waitbar((fr-st+1)/(en-st+1));

            msk_edge_now1 = double(imread(file_location1,fr));
            msk_edge_now2 = double(imread(file_location2,fr));
            msk_edge_now = zeros(size(msk_edge_now1));
            msk_edge_now(msk_edge_now1 + msk_edge_now2 > 0) = 1;
            im = imfill(msk_edge_now,'holes');
            im1 = imdilate(im,se);
            im1 = imfill(im1,'holes');
            im1 = imerode(im1,se);

            dis_edge_now = double(bwdist(~im1,'euclidean'));
            imdp = dis_edge_now - (dp-1);
            imdp(imdp<0) = 0; imdp(imdp>0) = 1;

            if fr == st
                ims = imdp;
                im_st = imdp;
            else
                ims = ims + imdp;
            end

        end
        close(wb);

        im_en = imdp;

        Xmax = size(ims,2);
        Ymax = size(ims,1);

        imin = ims;
        imin(imin<(T+1)) = 0;
        imin(imin>0) = 1;

        imou = ims;
        imou(imou>0) = 1;

        [Xin,Yin,~,~]=SmothBound(imin,1);
        [Xou,You,~,~]=SmothBound(imou,1);
        [Xst,Yst,~,~]=SmothBound(im_st,1);
        [Xen,Yen,~,~]=SmothBound(im_en,1);

        [Xin, Yin] = poly2ccw(Xin, Yin);
        [Xou, You] = poly2ccw(Xou, You);
        [Xst, Yst] = poly2ccw(Xst, Yst);
        [Xen, Yen] = poly2ccw(Xen, Yen);

        Xin = Xin + 0.0001*randn(size(Xin));
        Yin = Yin + 0.0001*randn(size(Yin));
        Xou = Xou + 0.0001*randn(size(Xou));
        You = You + 0.0001*randn(size(You));
        Xst = Xst + 0.0001*randn(size(Xst));
        Yst = Yst + 0.0001*randn(size(Yst));
        Xen = Xen + 0.0001*randn(size(Xen));
        Yen = Yen + 0.0001*randn(size(Yen));

        %% Finding middle loop 1

        X1 = Xin; Y1 = Yin;
        X2 = Xou; Y2 = You;
        LL = 1; cnt1 = 0;
        while LL > 0.01 && cnt1 < maxitr
            if cnt1 == 0 
                [Xm1,Ym1,Lm1] = mid_contour3a(X1,Y1,X2,Y2,['time block ' num2str(uu) ', mid-contour iteration ' num2str(cnt1+1) 'a']);
                Xm1o = Xm1; Ym1o = Ym1; Lm1o = Lm1;
                [Xm2,Ym2,Lm2] = mid_contour3a(X2,Y2,X1,Y1,['time block ' num2str(uu) ', mid-contour iteration ' num2str(cnt1+1) 'b']);
            else
                [Xm1,Ym1,Lm1] = mid_contour(X1,Y1,X2,Y2);
                [Xm2,Ym2,Lm2] = mid_contour(X2,Y2,X1,Y1);
            end

            X1 = Xm1; Y1 = Ym1;  
            X2 = Xm2; Y2 = Ym2; 
            LL = abs(Lm2-Lm1);
            cnt1 = cnt1 + 1;

        end

        XmdA1 = X1; YmdA1 = Y1; NmdA1 = length(XmdA1);
        %XmdB1 = X2; YmdB1 = Y2; NmdB1 = length(XmdB1);

        %% Finding middle loop 2

        X1 = Xin; Y1 = Yin;
        X2 = Xou; Y2 = You;
        LL = 1; cnt2 = 0;
        while LL > 0.01 && cnt2 < maxitr
            if cnt2 == 0  
                Xm1 = Xm1o; Ym1 = Ym1o; Lm1 = Lm1o;
                [Xm2,Ym2,Lm2] = mid_contour(X2,Y2,X1,Y1);
            else
                [Xm1,Ym1,Lm1] = mid_contour(X1,Y1,X2,Y2);
                [Xm2,Ym2,Lm2] = mid_contour(X2,Y2,X1,Y1);
            end

            X1 = Xm1; Y1 = Ym1;  
            X2 = Xm2; Y2 = Ym2; 
            LL = abs(Lm2-Lm1);
            cnt2 = cnt2 + 1;

        end

        XmdA2 = X1; YmdA2 = Y1; NmdA2 = length(XmdA2);
        %XmdB2 = X2; YmdB2 = Y2; NmdB2 = length(XmdB2);

        if cnt1<cnt2
            Xmd = XmdA1; Ymd = YmdA1; Nmd = NmdA1;
            convg(uu) = cnt1;
        else
            Xmd = XmdA2; Ymd = YmdA2; Nmd = NmdA2;
            convg(uu) = cnt2;
        end


        %% Finding normals

        Q1 = zeros(Nmd,2);
        Q2 = zeros(Nmd,2);
        Xs = zeros(4,Nmd); Ys = zeros(4,Nmd); Ds = zeros(4,Nmd);

        St = zeros(1,Nmd);

        wb = waitbar(0,['time block ' num2str(uu) ', normals']);
        for p = 1:Nmd

            waitbar(p/Nmd);

            if p ==1
                x1 = Xmd(Nmd); y1 = Ymd(Nmd);
                x2 = Xmd(p);   y2 = Ymd(p);
                x3 = Xmd(p+1); y3 = Ymd(p+1);
            elseif p == Nmd
                x1 = Xmd(p-1); y1 = Ymd(p-1);
                x2 = Xmd(p);   y2 = Ymd(p);
                x3 = Xmd(1);   y3 = Ymd(1);
            else
                x1 = Xmd(p-1); y1 = Ymd(p-1);
                x2 = Xmd(p);   y2 = Ymd(p);
                x3 = Xmd(p+1); y3 = Ymd(p+1);
            end

            if y1 == y3
                q1 = [1, y1]; q2 = [Xmax, y1];
            elseif x1 == x3
                q1 = [x1, 1]; q2 = [x1, Ymax];
            elseif abs((x1 - x3)/(y1 - y3)) < Ymax/Xmax  
                q1 = [1, y2 + (x1 - x3)*(x2 - 1)/(y1 - y3)];
                q2 = [Xmax, y2 + (x1 - x3)*(x2 - Xmax)/(y1 - y3)];
            else
                q1 = [x2 + (y1 - y3)*(y2 - 1)/(x1 - x3),1];
                q2 = [x2 + (y1 - y3)*(y2 - Ymax)/(x1 - x3), Ymax];
            end
            Q1(p,:) = q1;
            Q2(p,:) = q2;

            for w = 1:4
                Xcr = zeros(1,100);
                Ycr = zeros(1,100);
                Dcr = zeros(1,100);

                if w == 1
                    Xtm = Xin; Ytm = Yin;
                elseif w == 2
                    Xtm = Xou; Ytm = You;
                elseif w == 3
                    Xtm = Xst; Ytm = Yst;
                elseif w == 4
                    Xtm = Xen; Ytm = Yen;
                end    

                cnt = 0;
                r1 = [Xtm(end),Ytm(end)];
                r2 = [Xtm(1),Ytm(1)];
                [xi,yi,check] = intercheck(q1,q2,r1,r2);
                if check
                    cnt = cnt + 1; 
                    Xcr(cnt) = xi;
                    Ycr(cnt) = yi;
                    Dcr(cnt) = sqrt((x2-xi)^2+(y2-yi)^2);
                end
                for i = 2:length(Xtm)
                    r1 = [Xtm(i-1),Ytm(i-1)];
                    r2 = [Xtm(i),Ytm(i)];
                    [xi,yi,check] = intercheck(q1,q2,r1,r2);
                    if check
                        cnt = cnt + 1; 
                        Xcr(cnt) = xi;
                        Ycr(cnt) = yi;
                        Dcr(cnt) = sqrt((x2-xi)^2+(y2-yi)^2);
                    end
                end
                Xcr = Xcr(1:cnt);
                Ycr = Ycr(1:cnt);
                Dcr = Dcr(1:cnt);
                if isempty(Dcr)
                    St(p) = 1;
                    Ds(w,p) = 0;
                    Xs(w,p) =Xmd(p);
                    Ys(w,p) =Ymd(p);
                else
                    [Ds(w,p),j] = min(Dcr);
                    Xs(w,p) = Xcr(j); 
                    Ys(w,p) = Ycr(j);
                end
            end
            if ~St(p)
                x13 = x3 - x1;    y13 = y3 - y1;
                xio = Xs(2,p) - Xs(1,p);  yio = Ys(2,p) - Ys(1,p);
                if sign(y13*xio - x13*yio) == -1
                    St(p) = 1;
                    Ds(:,p) = 0;
                    Xs(:,p) =Xmd(p);
                    Ys(:,p) =Ymd(p);
                end
                    
            end

        end
        close(wb);
        
        cellXin{uu} = Xin;
        cellYin{uu} = Yin;
        cellXou{uu} = Xou;
        cellYou{uu} = You;
        cellQ1{uu} = Q1;
        cellQ2{uu} = Q2;
        cellXs{uu} = Xs;
        cellYs{uu} = Ys;
        cellDs{uu} = Ds;
        cellXmd{uu} = Xmd;
        cellYmd{uu} = Ymd;
        noINTER{uu} = St;
        N_noINT(uu) = sum(St);

        fig = figure('Position', [50 50 800 800]);
        hold on;
        axis ij; axis off; axis equal;
        %plot([Xmdw Xmdw(1)],[Ymdw Ymdw(1)],'c.-','LineWidth',2);
        plot([Xmd Xmd(1)],[Ymd Ymd(1)],'k.-','LineWidth',2);
        for p = 1:size(Xs,2)
            plot([Xs(1,p) Xs(2,p)],[Ys(1,p) Ys(2,p)],'k-');
        end
        plot([Xin(1,:) Xin(1,1)],[Yin(1,:) Yin(1,1)],'r.-');
        plot([Xou(1,:) Xou(1,1)],[You(1,:) You(1,1)],'b.-');
        %plot([Xs(3,:) Xs(3,1)],[Ys(3,:) Ys(3,1)],'ro');
        %plot([Xs(4,:) Xs(4,1)],[Ys(4,:) Ys(4,1)],'bo');
        
        plot(Xmd(St==1),Ymd(St==1),'go');
        
        drawnow;
        saveas(fig,fullfile(save_fold, strcat('block_', num2str(uu), '_dp', num2str(dp), '.fig')));
        saveas(fig,fullfile(save_fold, strcat('block_', num2str(uu), '_dp', num2str(dp), '.png')));
        close(fig);
        

        
    end
    
    save(fullfile(save_fold, strcat('normals_dp', num2str(dp), '.mat')), ...
        'cellXin','cellYin','cellXou','cellYou','cellQ1','cellQ2','cellXs','cellYs','cellDs','cellXmd','cellYmd','noINTER','convg','maxitr');
    %disp('total processing time (min):');
    elapsedTime = toc/60;
    %disp('number of non-interactions');
    %disp(N_noINT);
    %disp('number of steps before convergence');
    %disp(convg);
    