clear all
clc 

%%%%%%%%%%%%    EXAMPLE 1    %%%%%%%%%%%%
Points=[7.03 5.99;6.95 5.45;6.77 5.03;6.4 4.6;5.91 4.03;5.43 3.56;4.93 2.94;4.67 2.6;4.38 2.2;4.04 1.67;3.76 1.22;3.76 1.97;3.76 2.78;3.76 3.56;3.76 4.34;3.76 4.91;3.76 5.47;3.8 5.98;4.07 6.4;4.53 6.75; 5.07 6.85;5.05 6.84;5.89 6.83;6.41 6.8;6.92 6.58];
minV=[1 0 0 0 0 0 0 0 0];
maxV=[4 9 9 9 9 9 9 4 4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%    EXAMPLE 2    %%%%%%%%%%%%
% Points=[20 20;20 25;20 30;20 35;20 40;20 45];
% minV=[5 0 0 0 0 3.5  3.5  0  0];
% maxV=[55 60 60 60 60 38.5 38.5 60 60];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%    EXAMPLE 3    %%%%%%%%%%%%
% Points=[2.6 1.6;2.3 1.6;2 1.6;1.7 1.6;1.4 1.6;1 1.3;2 0.7;3 1.3];
% minV=[0 0 0 0 0 0 0 0 0];
% maxV=[5 4 4 4 4 3.5 3.5 4 4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
n=0;
ph(1:9)=0.5;
   d2=maxV(1)*rand()+minV(1);
   lx=((-1)^round(rand()))*maxV(2)*rand()+minV(2);
   ly=((-1)^round(rand()))*maxV(3)*rand()+minV(3);
   rx=((-1)^round(rand()))*maxV(4)*rand()+minV(4);
   ry=((-1)^round(rand()))*maxV(5)*rand()+minV(5);
   fx=((-1)^round(rand()))*maxV(6)*rand()+minV(6);
   fy=((-1)^round(rand()))*maxV(7)*rand()+minV(7);
   x0=((-1)^round(rand()))*maxV(8)*rand()+minV(8);
   y0=((-1)^round(rand()))*maxV(9)*rand()+minV(9);
   V=[d2 lx ly rx ry fx fy x0 y0];

old=zeros(1,9);
iterations=2000;
res=1;
exploration(1,1:9)=1;
while (i<iterations)
n=0;
sensitivity=sensitivityanalysis(V,Points);
    while (n<=9)
        if(n~=0)
            if(exploration(n)==0)
                   internaliteration=floor(abs(sensitivity(n)/10))+1;
            else
                internaliteration=1;
            end
        else
            internaliteration=1;
        end
        
        for j=1:internaliteration
            newtotalerror=0;
            deltaph=zeros(1,9);
            old=V;
%% changing parameter value            
        if(n ~= 0)
         if(exploration(n)==1)  % if exploration 
         if((ph(n))/(sum(ph))<rand())
         old(n)=V(n);
         if(n==1)
             V(n)=maxV(n)*rand()+minV(n);
         else
             V(n)=((-1)^round(rand()))*maxV(n)*rand()+minV(n);
         end
         end
         end
         if(exploration(n)==0) % if exploitation
         old(n)=V(n);
         V(n)=V(n)+((-1)^round(rand()))*rand()*(0.1*V(n)*(100-sensitivity(n))/100);
         
            if(V(n)>maxV(n)+minV(n))
                V(n)=maxV(n)+minV(n);
            end
            if(V(n)<minV(n))
                V(n)=minV(n);
            end
         end 
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% calculating error
totalerror=0;
error=zeros(360,size(Points,1));
minerror=zeros(1,size(Points,1));
l=sqrt(V(2)^2+V(3)^2);
r=sqrt(V(4)^2+V(5)^2);
d4=sqrt(V(6)^2+V(7)^2);
t4=atan2d(V(7),V(6));
 
% calculate d3
cx=V(2)-V(4);
cy=V(3)-V(5);
d3=sqrt(cx^2+cy^2);
 
% calculate t3
t3=atan2d(cy,cx);
A=atand(V(3)/V(2))-t3;
 
% Suppose t2 is 60 degrees, calculate d1 and t1
t2=60;
d1=sqrt(V(1)^2+d3^2+d4^2+2*V(1)*d3*cosd(t2-t3)-2*V(1)*d4*cosd(t2-t4)-2*d3*d4*cosd(t3-t4));
N1=d3*sind(t3)-d4*sind(t4)+V(1)*sind(t2);
N2=d3*cosd(t3)-d4*cosd(t4)+V(1)*cosd(t2);
t1=atan2d(N1,N2);
 
%check if it is not a rocker-rocker mechanism
D=[d1 V(1) d3 d4];
[small]=find(D==min(D));
[large]=find(D==max(D));
SL=D(small)+D(large);
PQ=sum(D)-D(small)-D(large);
if(SL<PQ && (D(small)==V(1)))
 for t2= 1:1:360

     %find t4
     u1=2*d1*d4*cosd(t1)-2*V(1)*d4*cosd(t2);
     u2=2*d1*d4*sind(t1)-2*V(1)*d4*sind(t2);
     u3=d1^2+V(1)^2-d3^2+d4^2-2*d1*V(1)*cosd(t1-t2);
     un=-u2-sqrt(u1^2+u2^2-u3^2);%changed sign before square root
     ud=u3-u1;
     t4=2*atan2d(un,ud);
     
     %find t3
     t3=atan2d(d1*sind(t1)+d4*sind(t4)-V(1)*sind(t2),d1*cosd(t1)+d4*cosd(t4)-V(1)*cosd(t2));
    
     %find the coordinates (xs,ys) of point S
     xs=V(1)*cosd(t2)+l*cosd(A+t3)+V(8);
     ys=V(1)*sind(t2)+l*sind(A+t3)+V(9);

    for m=1:size(Points,1)
    error(t2,m)=((xs-Points(m,1))^2)+((ys-Points(m,2))^2);
    end  
 end
for m=1:size(Points,1)
    minerror(m)=min(error(:,m));
    totalerror=totalerror+minerror(m);
end

newtotalerror=totalerror;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% comparing error value with the previous one and adding pheromones
if(n==0 && i==0)
     oldtotalerror=newtotalerror;
end
if(n~=0)
      if( newtotalerror<oldtotalerror)
      deltaph(n)=(oldtotalerror-newtotalerror)/oldtotalerror;
      ph(n)=0.2*ph(n)+deltaph(n);
      oldtotalerror=newtotalerror;
      else
          V(n)=old(n);
          ph(n)=0.2*ph(n);
      end
end
else
    if (n==0 && i==0)
        oldtotalerror=1000;
    end
    if(n==0  && i~=0)
        newtotalerror=oldtotalerror;
    end
    if(n~=0 && i~=0)
        newtotalerror=oldtotalerror;
        ph(n)=0.2*ph(n);
    end 
    if(n~=0)
    V(n)=old(n);
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% choosing next phase (exploitation or exploration)
   if (i==0)
      oldvalue=oldtotalerror; 
   end
    if (n~=0)
   if(i~=0)
   if(ph(n)<((0.2^7)*0.5))  % threshold
     exploration(n)=0;
     ph(n)=0.5;
   else
       oldvalue=oldtotalerror;
       exploration(n)=1;
   end
   end
    end
        end
    n=n+1;
    end  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   i=i+1;
   ITERATION=i
   V
   ERROR=oldtotalerror   
end