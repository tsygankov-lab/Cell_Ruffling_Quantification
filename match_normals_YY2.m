function [elapsedTime] = match_normals_YY2(fn1,dp,save_fold)

    tic

    load(fullfile(save_fold, strcat('normals_dp', num2str(dp), '.mat')));

    Nl = length(cellXs);
    L = zeros(1,Nl);

    for i = 1:Nl
        L(i) = length(cellXmd{i});
    end

    [Ntr,j] = max(L);

    IND_tr = zeros(Ntr,Nl);

    IND_tr(:,j) = (1:Ntr)';

    wb = waitbar(0,'matching progress...');
    for p = 1:Ntr

        waitbar(p/Ntr);

        for i = (j-1):-1:1
            Xs = cellXs{i+1};
            Ys = cellYs{i+1};
            St = noINTER{i+1};
            IND = IND_tr(p,i+1);
            st = St(IND);
            if IND == 0 || st == 1
                break;
            end
            X1 = Xs(3,IND);
            Y1 = Ys(3,IND);

            Xs = cellXs{i};
            Ys = cellYs{i};
            St = noINTER{i};
            X2 = Xs(4,:);
            Y2 = Ys(4,:);

            D = (X2 - X1).^2 + (Y2 - Y1).^2;
            D(St==1) = Inf;
            [mD,mI] = min(D);
            if mD>1
                break;
            end
            IND_tr(p,i) = mI;

        end

        for i = (j+1):Nl
            Xs = cellXs{i-1};
            Ys = cellYs{i-1};
            St = noINTER{i-1};
            IND = IND_tr(p,i-1);
            st = St(IND);
            if IND == 0 || st == 1
                break;
            end
            X1 = Xs(4,IND);
            Y1 = Ys(4,IND);

            Xs = cellXs{i};
            Ys = cellYs{i};
            St = noINTER{i};
            X2 = Xs(3,:);
            Y2 = Ys(3,:);

            D = (X2 - X1).^2 + (Y2 - Y1).^2;
            D(St==1) = Inf;
            [mD,mI] = min(D);
            if mD>1
                break;
            end
            IND_tr(p,i) = mI;

        end
    end
    close(wb);

    INDtr = IND_tr;
    cnt = 0;
    for p = 1:Ntr
       IND = IND_tr(p,:);
       f = find(IND==0);
       if isempty(f)
           cnt = cnt + 1;
           INDtr(cnt,:) = IND;
       end
    end
    Ntr = cnt;
    INDtr = INDtr(1:Ntr,:);
    
    wb = waitbar(0,'plotting result, please wait...');
    fig = figure('Position',get(0,'Screensize'),'visible','off');
    hold on;
    axis image; axis ij; axis off;
    for p = 1:Ntr
        waitbar(p/Ntr);
        
        for i = 1:Nl
            
            Xs = cellXs{i};
            Ys = cellYs{i};
            IND = INDtr(p,i);
            plot([Xs(1,IND),Xs(2,IND)],[Ys(1,IND),Ys(2,IND)],'k');
            
        end
    end
    close(wb);
    set(fig,'visible','on');
    
    saveas(fig,fullfile(save_fold, strcat('matched_dp', num2str(dp), '.fig')));
    saveas(fig,fullfile(save_fold, strcat('matched_dp', num2str(dp), '.png')));
    close(fig);
    
    save(fullfile(save_fold, strcat('matched_dp', num2str(dp), '.mat')),...
       'cellXin','cellYin','cellXou','cellYou','cellQ1','cellQ2','cellXs','cellYs','cellDs','cellXmd','cellYmd','noINTER','convg','maxitr','INDtr');
    
    
    elapsedTime = toc/60;
