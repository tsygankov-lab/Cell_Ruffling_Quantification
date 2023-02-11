function [elapsedTime] = quantify_along_normals_YY2(fn1,fn2,dp,save_fold)

    tic

    load(fullfile(save_fold, strcat('matched_dp', num2str(dp), '.mat')));

    file_location1 = fn1;
    file_location2 = fn2;

    
    info = imfinfo(file_location1);
    %N_frames = length(info);
    T = 20;
    cont_N = fix(length(info)/T);
    if mod(length(info),T) == 0
        cont_N = cont_N - 1;
    end
    N_frames = 1 + T*cont_N;

    se = strel('disk',1);

    Xcnt = cell(1,N_frames);
    Ycnt = cell(1,N_frames);

    wb = waitbar(0,'processing contours...');
    for fr = 1:N_frames

            waitbar(fr/N_frames);

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

            if fr == 1
                ims = imdp;
            else
                ims = ims + imdp;
            end

            [Xt,Yt,~,~]=SmothBound(imdp,1);
            Xcnt{fr} =  Xt + 0.0001*randn(size(Xt));
            Ycnt{fr} =  Yt + 0.0001*randn(size(Yt));


    end
    close(wb);

    N_tr = size(INDtr,1);
    x_tr = zeros(N_frames,N_tr);
    y_tr = zeros(N_frames,N_tr);
    v_tr = zeros(N_frames-1,N_tr);

    fig = figure('Position',get(0,'Screensize'),'visible','off');
    hold on;
    axis image; axis ij; axis off;
    wb = waitbar(0,'finding velocity values...');
    for p = 1:N_tr

        waitbar(p/N_tr);

        for uu = 1:cont_N

            st = 1+T*(uu-1);
            en = 1+T*uu;

            IND = INDtr(p,uu);
            Xs = cellXs{uu};
            Ys = cellYs{uu};
            Xmd = cellXmd{uu};
            Ymd = cellYmd{uu};

            Q1 = cellQ1{uu};
            Q2 = cellQ2{uu};

            q1 =Q1(IND,:); q2 =Q2(IND,:);

            for fr = st:en%1:N_frames
                X = Xcnt{fr}; Y = Ycnt{fr};

                Xcr = zeros(1,100);
                Ycr = zeros(1,100);
                Dcr = zeros(1,100);
                cnt = 0;

                r1 = [X(end),Y(end)];
                r2 = [X(1),Y(1)];
                [xi,yi,check] = intercheck(q1,q2,r1,r2);
                if check
                    cnt = cnt + 1; 
                    Xcr(cnt) = xi;
                    Ycr(cnt) = yi;
                    Dcr(cnt) = sqrt((Xmd(IND)-xi)^2+(Ymd(IND)-yi)^2);
                end
                for i = 2:length(X)
                    r1 = [X(i-1),Y(i-1)];
                    r2 = [X(i),Y(i)];
                    [xi,yi,check] = intercheck(q1,q2,r1,r2);
                    if check
                        cnt = cnt + 1; 
                        Xcr(cnt) = xi;
                        Ycr(cnt) = yi;
                        Dcr(cnt) = sqrt((Xmd(IND)-xi)^2+(Ymd(IND)-yi)^2);
                    end
                end
                Xcr = Xcr(1:cnt);
                Ycr = Ycr(1:cnt);
                Dcr = Dcr(1:cnt);
                [~,j] = min(Dcr);

                if isempty(j)
                    continue;
                end

                x_tr(fr,p) = Xcr(j); y_tr(fr,p) = Ycr(j);

                if fr>1
                    d = (x_tr(fr,p) - x_tr(fr-1,p)).^2+(y_tr(fr,p) - y_tr(fr-1,p)).^2;

                    if inpolygon(x_tr(fr-1,p),y_tr(fr-1,p),X,Y)
                        v_tr(fr-1,p) = sqrt(d);   
                    else
                        v_tr(fr-1,p) = -sqrt(d);
                    end
                end

            end

        end

        plot(x_tr(:,p),y_tr(:,p),'b-');

    end
    close(wb);
    set(fig,'visible','on');
    
    saveas(fig,fullfile(save_fold, strcat('tracked_dp', num2str(dp), '.fig')));
    saveas(fig,fullfile(save_fold, strcat('tracked_dp', num2str(dp), '.png')));
    close(fig);
    
    save(fullfile(save_fold, strcat('tracked_dp', num2str(dp), '.mat')),'x_tr','y_tr','v_tr','Xcnt','Ycnt');
    
    elapsedTime = toc/60;