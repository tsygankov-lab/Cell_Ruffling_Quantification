function [xi,yi,check] = intercheck(q1,q2,r1,r2)

    check = 0; xi = []; yi = [];
    
    eps = 1e-9;
    
    if q1(1)==r1(1) && q1(2)==r1(2) && q2(1)==r2(1) && q2(2)==r2(2) 
        return;
    elseif q1(1)==r2(1) && q1(2)==r2(2) && q2(1)==r1(1) && q2(2)==r1(2)
        return;
    end

    if q1(1)~=q2(1)

        kq = (q2(2)-q1(2))/(q2(1)-q1(1));
        if q1(1) == 0
            bq = q1(2);
        elseif q2(1) == 0
            bq = q2(2);
        else
            bq = (q2(1)*q1(2)-q1(1)*q2(2))/(q2(1)-q1(1));
        end

        if r1(1)~=r2(1)

            kr = (r2(2)-r1(2))/(r2(1)-r1(1));
            if r1(1) == 0
                br = r1(2);
            elseif r2(1) == 0
                br = r2(2);
            else
                br = (r2(1)*r1(2)-r1(1)*r2(2))/(r2(1)-r1(1));
            end
            
            if kq~=kr

                xi = - (bq - br)/(kq - kr);
                if kq == 0
                    yi = bq;
                elseif kr == 0
                    yi = br;
                else
                    yi = (kq*br - kr*bq)/(kq - kr);
                end

                %if xi>=min(q1(1),q2(1)) && xi<=max(q1(1),q2(1)) && yi>=min(q1(2),q2(2)) && yi<=max(q1(2),q2(2))
                %    if xi>=min(r1(1),r2(1)) && xi<=max(r1(1),r2(1)) && yi>=min(r1(2),r2(2)) && yi<=max(r1(2),r2(2))
                if (xi - min(q1(1),q2(1)))>-eps && (max(q1(1),q2(1))-xi)>-eps && (yi-min(q1(2),q2(2)))>-eps && (max(q1(2),q2(2))-yi)>-eps
                    if (xi - min(r1(1),r2(1)))>-eps && (max(r1(1),r2(1))-xi)>-eps && (yi-min(r1(2),r2(2)))>-eps && (max(r1(2),r2(2))-yi)>-eps
                        check = 1;
                    end
                end
            end
        else
            xi = r1(1);
            yi = kq*xi + bq;
            %if xi>=min(q1(1),q2(1)) && xi<=max(q1(1),q2(1)) && yi>=min(q1(2),q2(2)) && yi<=max(q1(2),q2(2))
            %    if xi>=min(r1(1),r2(1)) && xi<=max(r1(1),r2(1)) && yi>=min(r1(2),r2(2)) && yi<=max(r1(2),r2(2))
            if (xi - min(q1(1),q2(1)))>-eps && (max(q1(1),q2(1))-xi)>-eps && (yi-min(q1(2),q2(2)))>-eps && (max(q1(2),q2(2))-yi)>-eps
                if (xi - min(r1(1),r2(1)))>-eps && (max(r1(1),r2(1))-xi)>-eps && (yi-min(r1(2),r2(2)))>-eps && (max(r1(2),r2(2))-yi)>-eps
                    check = 1;
                end
            end
        end
    elseif r1(1)~=r2(1)

        kr = (r2(2)-r1(2))/(r2(1)-r1(1));
        if r1(1) == 0
            br = r1(2);
        elseif r2(1) == 0
            br = r2(2);
        else
            br = (r2(1)*r1(2)-r1(1)*r2(2))/(r2(1)-r1(1));
        end
        
        xi = q1(1);
        yi = kr*xi + br;

        %if xi>=min(q1(1),q2(1)) && xi<=max(q1(1),q2(1)) && yi>=min(q1(2),q2(2)) && yi<=max(q1(2),q2(2))
        %    if xi>=min(r1(1),r2(1)) && xi<=max(r1(1),r2(1)) && yi>=min(r1(2),r2(2)) && yi<=max(r1(2),r2(2))
        if (xi - min(q1(1),q2(1)))>-eps && (max(q1(1),q2(1))-xi)>-eps && (yi-min(q1(2),q2(2)))>-eps && (max(q1(2),q2(2))-yi)>-eps
            if (xi - min(r1(1),r2(1)))>-eps && (max(r1(1),r2(1))-xi)>-eps && (yi-min(r1(2),r2(2)))>-eps && (max(r1(2),r2(2))-yi)>-eps
                check = 1;
            end
        end   
    end

    %if check
    %    figure;
    %    hold on;
    %    plot([q1(1) q2(1)],[q1(2) q2(2)],'r','LineWidth',3);
    %    plot([r1(1) r2(1)],[r1(2) r2(2)],'b','LineWidth',3);
    %    if check
    %        plot(xi,yi,'ko','MarkerFaceColor','y');
    %    end
    %end

end



