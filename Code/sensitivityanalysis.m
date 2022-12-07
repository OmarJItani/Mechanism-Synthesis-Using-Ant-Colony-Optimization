function [sensitivity]=sensitivityanalysis(V,Points)
newV=V;
h=zeros(1,9);
for m=1:10
    totalerror=0;
    if(m>1)
        h(m-1)=0.05*V(m-1);
        newV(m-1)=V(m-1)+h(m-1);
    end
    d2=newV(1);
    l=sqrt(newV(2)^2+newV(3)^2);
    r=sqrt(newV(4)^2+newV(5)^2);
    d4=sqrt(newV(6)^2+newV(7)^2);
    t4=atan2d(newV(7),newV(6));

    % calculate d3
    cx=newV(2)-newV(4);
    cy=newV(3)-newV(5);
    d3=sqrt(cx^2+cy^2);

    % calculate t3 
    t3=atan2d(cy,cx);
    A=atand(newV(3)/newV(2))-t3;

    % Suppose t2 is 60 degrees, calculate d1 and t1
    t2=60;
    d1=sqrt(newV(1)^2+d3^2+d4^2+2*newV(1)*d3*cosd(t2-t3)-2*newV(1)*d4*cosd(t2-t4)-2*d3*d4*cosd(t3-t4));
    N1=d3*sind(t3)-d4*sind(t4)+newV(1)*sind(t2);
    N2=d3*cosd(t3)-d4*cosd(t4)+newV(1)*cosd(t2);
    t1=atan2d(N1,N2);
    
    %check if it is not a rocker-rocker mechanism
    D=[d1 newV(1) d3 d4];
    [small]=find(D==min(D));
    [large]=find(D==max(D));
    SL=D(small)+D(large);
    PQ=sum(D)-D(small)-D(large);
    if(SL<PQ && (D(small)==newV(1)))
        for t2= 1:360
     
         %find t4
         u1=2*d1*d4*cosd(t1)-2*newV(1)*d4*cosd(t2);
         u2=2*d1*d4*sind(t1)-2*newV(1)*d4*sind(t2);
         u3=d1^2+newV(1)^2-d3^2+d4^2-2*d1*newV(1)*cosd(t1-t2);
         un=-u2-sqrt(u1^2+u2^2-u3^2);
         ud=u3-u1;
         t4=2*atan2d(un,ud);  
         %find t3
         t3=atan2d(d1*sind(t1)+d4*sind(t4)-newV(1)*sind(t2),d1*cosd(t1)+d4*cosd(t4)-newV(1)*cosd(t2));  
         %find the coordinates (xs,ys) of point S
         xs=newV(1)*cosd(t2)+l*cosd(A+t3)+newV(8);
         ys=newV(1)*sind(t2)+l*sind(A+t3)+newV(9);
         
         for j=1:size(Points,1)
         error(t2,j)=((xs-Points(j,1))^2)+((ys-Points(j,2))^2);
         end  
        end
    for i=1:size(Points,1)
    minerror(i)=min(error(:,i));
    totalerror=totalerror+minerror(i);
    end
    else
        totalerror=1000;
    end
    sensitivityerror(m)=totalerror;
    if(m>1)
    newV(m-1)=newV(m-1)-h(m-1);
    end
end
for k=2:10
    sensitivity(k-1)=(sensitivityerror(k)-sensitivityerror(1))/h(k-1);
end
sensitivity=abs(sensitivity);
[p]=find(sensitivity==max(sensitivity));
sensitivity=abs(round((sensitivity*100)/sensitivity(p(1))));
end