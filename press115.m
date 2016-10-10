function data=press115(alph,bet,gam,D,n_owe,g);
% press115 calculates data from input parameters
% contains list of points from x-ray images


h=6.626*10^-34;
sp=2.998*10^8;
e=1.60217733*10^-19;
V=35*10^3;
pi=3.14159265358979;
lam_min=h*sp/(e*V);
lam_max=1.38059*10^(-10);

%LATTICE CONSTANTS

a=4.6662*10^(-10);
b=4.6662*10^(-10);
c=7.5168*10^(-10);
% D=3.44*10^(-2);
% %At zero pressure
% %a=4.6662*10^(-10);
% %b=4.6662*10^(-10);
% %c=7.5168*10^(-10);
% %The Indium atoms are located at sites: (0, .5, .3053) (like pressure
% %induced CeRhIn5; Yoshinori Haga, et al. Phys. Rev. B, 63, 060503, 2001.
% 
% %CHOOSE INCIDENT BEAM DIRECTION
% 
% % When analyzing the pressure effects, graphing with the c-axis into the
% % plane of the image, the "a" dimension runs vertically and the "b"
% % dimension runs horizontally.
% % max n_o||x,y is 26; n_o || z is 42
% dink=1;
% while dink==1
%     n_owe=-1*input('input the incident beam direction: only either [1,0,0] , [0,1,0] , [0,0,1] allowed:  ');
%     
% %     [1 0 0] has c-axis running horizontally, [0 1 0] c-axis runs vertically
%     
%     if sum(n_owe)~=-1 | (n_owe*n_owe.')~=1
%         display('try again');
%     elseif n_owe(1)==-1 & n_owe(2)==0 & n_owe(3)==0
%         dink=2;
%     elseif n_owe(1)==0 & n_owe(2)==-1 & n_owe(3)==0
%         dink=2;
%     elseif n_owe(1)==0 & n_owe(2)==0 & n_owe(3)==-1
%         dink=2;
%     else
%         display('try again');
%     end
% end
% 
% %INPUT ANGLES
% 
% alph=(pi/180)*input('input angle alpha in degrees, see code for angle definition:  ');
% bet=(pi/180)*input('input angle beta in degrees:  ');
% gam=(pi/180)*input('input angle gamma in degrees:  ');
% %The angle alpha is the rotation about the incident beam direction such
% %that, as viewed from that direction.  Beta is the rotation about the crystal axis aligned with the "y" axis on the image when alpha =0.  
% %Gamma is the rotation about the crystal axis aligned with the "x" axis on the image when alpha =0.  Of course the order of
% %operation of the rotation matrices affects the resultant crystal
% %orientation.  They are set to operate in a particular order:  1st=Gamma,
% %2nd=Beta, 3rd=Alpha.  However, one nonzero angle at a time is best and the
% %rotation becomes the identity matrix when an angle is zero.
% 
% %INPUT MAXIMUM HKL VALUE
% 
% g=input('type the abs(maximum) hkl value you wish to be evaluated;optimum=[0 0 1]-->30,[0 1 0] | [1 0 0]-->20:  ');

%ROTATION MATRICES DEFINED

n_o=n_owe.';
if alph~=0 | bet~=0 | gam~=0
    if n_o(3)~=0
        Ro=[cos(alph) -sin(alph) 0;sin(alph) cos(alph) 0; 0 0 1];
        Ry=[1 0 0;0 cos(bet) -sin(bet);0 sin(bet) cos(bet)];
        Rx=[cos(gam) 0 -sin(gam);0 1 0;sin(gam) 0 cos(gam)];
    elseif n_o(2)~=0
        Ry=[cos(bet) -sin(bet) 0;sin(bet) cos(bet) 0; 0 0 1];
        Rx=[1 0 0;0 cos(gam) -sin(gam);0 sin(gam) cos(gam)];
        Ro=[cos(alph) 0 -sin(alph);0 1 0;sin(alph) 0 cos(alph)];
    elseif n_o(1)~=0
        Rx=[cos(gam) -sin(gam) 0;sin(gam) cos(gam) 0; 0 0 1];
        Ro=[1 0 0;0 cos(alph) -sin(alph);0 sin(alph) cos(alph)];
        Ry=[cos(bet) 0 -sin(bet);0 1 0;sin(bet) 0 cos(bet)];
    end
end

%INITIALIZE VARIOUS QUANTITIES

cnt_bad=0;
cnt=0;
over=0;
under=0;

if g>=12 & g<16
    ourdata=zeros(4500,10);
    bad_dates=zeros(9500,8);
    data=zeros(750,10);
elseif g>=16 & g<20
    if n_o(1)~=0 || n_o(2)~=0
        ourdata=zeros(9000,10);
        bad_dates=zeros(26000,8);
        data=zeros(2500,10);
    elseif n_o(3)~=0
        ourdata=zeros(9000,10);
        bad_dates=zeros(26000,8);
        data=zeros(600,10);
    end
elseif g>=20 & g<24
    ourdata=zeros(12500,10);
    bad_dates=zeros(55000,8);
    data=zeros(5000,10);
elseif g>=24
    ourdata=zeros(14000,10);
    bad_dates=zeros(102000,8);
    data=zeros(6500,10);
end

%COMPUTE VIABLE POLES

if alph~=0 | bet~=0 | gam~=0
    for m3=-g:g;
        for m2=-g:g
            for m1 = -g:g;
                
                if m3==0 & m2==0 & m1==0
                    
                else
                    enn=Ro*Ry*Rx*[m1/a;m2/b;m3/c];
                    lam0=-2*(enn(1)*n_o(1)+enn(2)*n_o(2)+enn(3)*n_o(3))/(enn(1)^2+enn(2)^2+enn(3)^2);
                    if lam0>=lam_min && lam0<=lam_max
                        cnt=cnt+1;
                        theta=(180/pi)*asin(.5*lam0*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        ourdata(cnt,:)=[m1 m2 m3 real(n) real(theta) real(lam0) real(our_vbl) real(mag_n)];
                    elseif lam0< lam_min
                        under=under+1;
                        cnt_bad=cnt_bad+1;
                        theta=(180/pi)*asin(.5*lam_min*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        bad_dates(cnt_bad,:)=[m1 m2 m3 real(n) real(mag_n) real(n_dev)];
                    elseif lam0>lam_max
                        over=over+1;
                        cnt_bad=cnt_bad+1;
                        theta=(180/pi)*asin(.5*lam_max*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        bad_dates(cnt_bad,:)=[m1 m2 m3 real(n) real(mag_n) real(n_dev)];
                    end
                end
                        

                

                
            end
        end
    end
else
    
    % COMPUTER POLES FOR ZERO SET OF ANGLES (QUICKER)
    
    for m3=-g:g;
        for m2=-g:g
            for m1 = -g:g;
                
                if m3==0 & m2==0 & m1==0
                    
                else
                    enn=[m1/a;m2/b;m3/c];
                    lam0=-2*(enn(1)*n_o(1)+enn(2)*n_o(2)+enn(3)*n_o(3))/(enn(1)^2+enn(2)^2+enn(3)^2);
                    if lam0>=lam_min && lam0<=lam_max
                        cnt=cnt+1;
                        theta=(180/pi)*asin(.5*lam0*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        ourdata(cnt,:)=[m1 m2 m3 real(n) real(theta) real(lam0) real(our_vbl) real(mag_n)];
                    elseif lam0< lam_min
                        under=under+1;
                        cnt_bad=cnt_bad+1;
                        theta=(180/pi)*asin(.5*lam_min*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        bad_dates(cnt_bad,:)=[m1 m2 m3 real(n) real(mag_n) real(n_dev)];
                    elseif lam0>lam_max
                        over=over+1;
                        cnt_bad=cnt_bad+1;
                        theta=(180/pi)*asin(.5*lam_max*sqrt(enn.'*enn));
                        our_vbl=.5*sqrt(enn.'*enn);
                        n=[lam0*enn(1)+n_o(1),lam0*enn(2)+n_o(2),lam0*enn(3)+n_o(3)];
                        mag_n=(n(1)^2+n(2)^2+n(3)^2)^(.5);
                        n_dev=mag_n-1;
                        bad_dates(cnt_bad,:)=[m1 m2 m3 real(n) real(mag_n) real(n_dev)];
                    end
                end
            end
            
        end
    end
end

%SORT DATA

ord_data=sortrows(ourdata,9);   %sorted by our_vbl,
ord_bad=sortrows(bad_dates,-7);

%DETERMINE (X,Y) POSITION (LIMITED BY ANGLE)

if n_o(1)~=0
    blip=0;
    for k=1:size(ord_data,1)
        if ord_data(k,7)>=57.36
            blip=blip+1;
            r(blip,:)=100*D*[ord_data(k,6)./abs(ord_data(k,4)),ord_data(k,5)./abs(ord_data(k,4))];
            data(blip,:)=ord_data(k,:);
        end
    end
elseif n_o(2)~=0
    blip=0;
    for k=1:size(ord_data,1)
        if ord_data(k,7)>=57.36
            blip=blip+1;
            r(blip,:)=100*D*[ord_data(k,4)./abs(ord_data(k,5)),ord_data(k,6)./abs(ord_data(k,5))];
            data(blip,:)=ord_data(k,:);
        end
    end
elseif n_o(3)~=0
    blip=0;
    for k=1:size(ord_data,1)
        if ord_data(k,7)>=57.36
            blip=blip+1;
            r(blip,:)=100*D*[ord_data(k,5)./abs(ord_data(k,6)),ord_data(k,4)./abs(ord_data(k,6))];
            data(blip,:)=ord_data(k,:);
        end
    end
end
data(:,12:13)=real(r);

%DETERMINE ATOMIC FORM FACTORS

vbl_dsct(:,1)=(10^10)*[0 .01 .02 .03 .04 .05 .06 .07 .08 .09 .1 .11 .12 .13 .14 .15 .16 .17 .18 .19 .20 .22 .24 .25 .26 .28 .30 .32 .34 .35 .36 .38 .40 .42 .44 .45 .46 .48 .50 .55 .60 .65 .70 .80 .90 1.00 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80 1.90 2.00].';
f_Ce=[58 57.928 57.715 57.375 56.924 56.385 55.779 55.127 54.446 53.750 53.047 52.345 51.646 50.952 50.263 49.579 48.901 48.227 47.557 46.892 46.233 44.933 43.663 43.042 42.432 41.244 40.104 39.014 37.975 37.474 36.985 36.040 35.139 34.277 33.451 33.051 32.658 31.893 31.154 29.409 27.791 26.289 24.901 22.469 20.481 18.881 17.583 16.491 15.526 14.633 13.776 12.939 12.123 11.333 10.576 9.868].';
f_In=[49 48.957 48.828 48.618 48.332 47.980 47.570 47.112 46.614 46.086 45.534 44.964 44.383 43.793 43.199 42.603 42.006 41.410 40.817 40.226 39.639 38.478 37.337 36.774 36.218 35.125 34.059 33.025 32.025 31.538 31.060 30.134 29.247 28.401 27.596 27.209 26.832 26.108 25.425 23.881 22.552 21.405 20.408 18.736 17.329 16.053 14.840 13.670 12.548 11.492 10.518 9.639 8.860 8.184 7.603 7.110].';
f_Ir=[77 76.951 76.806 76.567 76.240 75.832 75.352 74.806 74.206 73.558 72.872 72.156 71.416 70.658 69.887 69.108 68.324 67.538 66.752 65.969 65.189 63.645 62.127 61.380 60.641 59.189 57.773 56.395 55.056 54.401 53.756 52.496 51.274 50.091 48.946 48.387 47.837 46.765 45.726 43.269 40.994 38.878 36.901 33.305 30.125 27.323 24.882 22.789 21.019 19.541 18.312 17.287 16.422 15.678 15.018 14.418].';

f1=zeros(size(data,1),1);
f2=zeros(size(data,1),1);
f3=zeros(size(data,1),1);
for j=1:size(data,1)
    vbl_min(1)=vbl_dsct(1);
    vbl_max(1)=vbl_dsct(end);
    
    if data(j,9)>vbl_min(1) && data(j,9)<vbl_max(1)
       bng=0;
       plc1(1)=1;
       plc2(1)=size(vbl_dsct,1);
       while plc2(bng+1)-plc1(bng+1)>1
            bng=bng+1;
            plc2(bng+1)=ceil(.5*(plc1(bng)+plc2(bng)));
            vbl_max(bng+1)=vbl_dsct(plc2(bng+1));
            if data(j,9)<vbl_max(bng+1)
                plc1(bng+1)=floor(.5*(plc1(bng)+plc2(bng+1)));
                vbl_min(bng+1)=vbl_dsct(plc1(bng+1));
                if data(j,9)<vbl_min(bng+1)
                    plc2(bng+1)=plc1(bng+1);
                    plc1(bng+1)=plc1(bng);
                    vbl_max(bng+1)=vbl_dsct(plc2(bng+1));
                    vbl_min(bng+1)=vbl_dsct(plc1(bng+1));
                    continue
                elseif data(j,9)==vbl_min(bng+1)
                    plc2(bng+1)=plc1(bng+1);
                    vbl_max(bng+1)=vbl_min(bng+1)
                    break
                end
            elseif data(j,9)>vbl_max(bng+1)
                plc1(bng+1)=plc2(bng+1);
                plc2(bng+1)=plc2(bng);
                vbl_min(bng+1)=vbl_dsct(plc1(bng+1));
                vbl_max(bng+1)=vbl_dsct(plc2(bng+1));
                continue
            elseif data(j,9)==vbl_max(bng+1)
                plc1(bng+1)=plc2(bng+1);
                vbl_min(bng+1)=vbl_max(bng+1);
                break
            end
       end
       if plc2(bng+1)==plc1(bng+1)
           f1(j)=f_Ce(plc2(bng+1));
           f2(j)=f_Ir(plc2(bng+1));
           f3(j)=f_In(plc2(bng+1));
           continue
       elseif plc2(bng+1)-plc1(bng+1)==1
           f1(j)=f_Ce(plc2(bng+1))+((f_Ce(plc2(bng+1))-f_Ce(plc1(bng+1)))/(vbl_dsct(plc2(bng+1))-vbl_dsct(plc1(bng+1))))*(data(j,9)-vbl_dsct(plc2(bng+1)));
           f2(j)=f_Ir(plc2(bng+1))+((f_Ir(plc2(bng+1))-f_Ir(plc1(bng+1)))/(vbl_dsct(plc2(bng+1))-vbl_dsct(plc1(bng+1))))*(data(j,9)-vbl_dsct(plc2(bng+1)));
           f3(j)=f_In(plc2(bng+1))+((f_In(plc2(bng+1))-f_In(plc1(bng+1)))/(vbl_dsct(plc2(bng+1))-vbl_dsct(plc1(bng+1))))*(data(j,9)-vbl_dsct(plc2(bng+1)));
       elseif plc2(bng+1)-plc1(bng+1)<1
           display('error2')
           return
       end
   elseif data(j,9)<vbl_min(1)
       display('error1')
       return
   elseif data(j,9)>vbl_max(1)
       f1(j)=f_Ce(end)+(f_Ce(end)-f_Ce(end-1))/(vbl_dsct(end)-vbl_dsct(end-1))*(data(j,9)-vbl_dsct(end,1));
       f2(j)=f_Ir(end)+(f_Ir(end)-f_Ir(end-1))/(vbl_dsct(end)-vbl_dsct(end-1))*(data(j,9)-vbl_dsct(end,1));
       f3(j)=f_In(end)+(f_In(end)-f_In(end-1))/(vbl_dsct(end)-vbl_dsct(end-1))*(data(j,9)-vbl_dsct(end,1));
   end
end
evod=data(:,3)/2;
evener=data(:,2)/4;
Trm3=zeros(size(data,1),1);
for t=1:size(data,1)
    if data(t,1)==0 & data(t,2)==0 & data(t,3)==0
        display('zuh?: [m1 m2 m3]=[0 0 0]')
        [data(t,1) data(t,2) data(t,3)]
        display('yep')
    elseif data(t,1)~=0 | data(t,2)~=0 | data(t,3)~=0
        ms2a(t,:)=[data(t,1) data(t,2) data(t,3)];
        if mod(evod(t),2)~=0
            Trm3(t)=(-1)^(data(t,1)+data(t,2));
        elseif mod(evener(t),2)==0
            Trm3(t)=(-1)^(data(t,1)+data(t,2))+2*((-1)^(data(t,1)+(-1)^(data(t,2))));
        elseif mod(evener(t),2)~=0
            Trm3(t)=(-1)^(data(t,1)+data(t,2))-2*((-1)^(data(t,1)+(-1)^(data(t,2))));
        end
    end
end

%DEFINE INTENSITY

Sk=f1+f2.*((-1).^(data(:,3)))+f3.*Trm3;
Ints=conj(Sk).*Sk;
data(:,11)=Ints(:);
check=sortrows(ord_data,-7);
check_bad=sortrows(ord_bad,3);
check1=sortrows(data,-12);
check2=sortrows(data,-13);
data1=sortrows(data,-3);
data2=sortrows(data1,-2);
data3=sortrows(data2,-1);
data=data3;

figure(1)

% xy1a: [angles]=[0,-.5,-1] [1 0 0]
% xy1a=[.036,3.762
%         0,2.124
%         -1.876,2.914
%         1.896,2.856
%         -2.118 ,-.070
%         2.070,-.142
%         -1.454,2.250
%         1.492,2.242
%         -1.670,1.234
%         1.684,1.186
%         -2.226,1.684
%         2.246,1.650
%         -2.736,1.372
%         2.766,1.330
%         -3.438,1.040
%         3.452,.926
%         -3.714,-.054
%         3.660,-.180
%         -1.746,-1.476
%         1.720,-1.538
%         -2.884,-1.696
%         2.818,-1.680
%         -2.364,-1.950
%         2.286,-2.026
%         -0.038,-2.492
%         -0.088,-4.432
%         -2.536,-2.774
%         2.398,-2.796
%         -3.440,2.716[
%         3.440,2.680
%         -3.590,-1.240
%         3.518,-1.320
%         -1.530,-3.628
%         1.404,-3.638];
%     plot(r(:,1),r(:,2),'.b',xy1a(:,1),xy1a(:,2),'*r')
%     axis([-4.4 4.4 -5.7 5.7])

% xy1b: [angles]=[0,-1.8,4] [0 0 1]
% xy1b=[.230,2.280
%         -2.656,-.486
%         3.430,-.520
%         .236,-3.866
%         -2.880,2.566
%         3.720,2.884
%         .214,-2.898
%         1.690,-1.948
%         1.488,-1.712
%         -1.216,-1.860
%         -1.196,2.318
%         1.668,2.396
%         -2.636,.948
%         3.382,1.084
%         -2.616,.468
%         3.358,.544
%         1.470,-4.092
%         2.258,-4.29
%         -.876,-3.92
%         -1.464,-4.026];
% plot(r(:,1),r(:,2),'.b',xy1b(:,1),xy1b(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

% xy2: [angles]=[-1.9,-2.4,.75] [0 0 1]
% xy2=[-2.532,-.134
%     .220,2.782
%     -2.602,.786
%     .410,-3.134
%     -2.588,-1.148
%     -2.658,-1.650
%     1.218,3.018
%     1.814,3.184
%     -.828,2.848
%     -1.302,2.880
%     -.668,-3.172
%     1.488,-3.256
%     3.518,.042
%     3.536,1.130
%     3.670,-1.088];
%  plot(r(:,1),r(:,2),'.b',xy2(:,1),xy2(:,2),'*r')
%  axis([-4.4 4.4 -5.7 5.7])

% %xy3b:  [angles]=[-1.5,-4,-5]   [1 0 0]
% xy3b=[-2.76,-.78
%     1.60,-.59
%     -.48,1.50
%     1.26,3.98];
% plot(r(:,1),r(:,2),'.b',xy3b(:,1),xy3b(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

%xy3c:  [angles]=[0 2 4]   [1 0 0]
% xy3c=[-1.87,.44
%     2.38,.40
%     .21,2.17
%     3.34,2.15];
% plot(r(:,1),r(:,2),'.b',xy3c(:,1),xy3c(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

% %xy3d:  [angles]=[   [1 0 0]
%xy3d=[-2.06,-.53
 %   2.21,-.72
  %  .11,1.52
%    3.07,-2.59];
%plot(r(:,1),r(:,2),'.b',xy3d(:,1),xy3d(:,2),'*r')
%axis([-4.4 4.4 -5.7 5.7])

%xy4:[angles]=[.25,1.6,1.25] [0 0 1]
% xy4=[-3.34,-.14
%     -4.07,-.17
%     2.68,-.16
%     3.29,-.19
%     -.26,-3.226
%     -.25,3.308
%     .51,3.17
%     -.9,3.2
%     -3.36,.86
%     2.7,.69
%     2.8,-1.15
%     -3.4,-1.26
%     -1.31,-3.4
%     .834,-3.270
%     -.2,-2.39
%     -.99,-2.55
%     .6,-2.45
%     -.180,2.684
%     .806,2.716
%     -1.168,2.802
%     -1.740,2.898
%     1.298,2.774];
% xy=xy4;
% 
%xy5: [angles]=[6.4,-2.25,1.25] [0 0 1]
% xy5=[.538,2.752
%     1.604,2.772
%     -.432,2.860
%     .666,3.308
%     -2.540,.162
%     -3.080,.240
%     -2.446,1.124
%     -2.696,-.806
%     -3.114,-.884
%     -2.822,1.252
%     3.464,-.514
%     3.654,.606
%     3.466,-1.560
%     -.116,-3.886
%     -.034,-3.184
%     1.042,-3.448
%     -1.072,-3.070];
% xy=xy5;
% 
%xy5_1: D=3.4875cm [angles]=[0 5.6 -.5] [0 1 0]
% xy5_1a=2.54*[-.299,.907
%     .303,.681
%     .679,.571
%     .989,.509
%     .660,-.027
%     1.218,-.029
%     .493,.290
%     .672,.368
%     .317,1.127
%     .446,.896
%     1.063,.319
%     .043,1.020
%     .215,.883
%     .723,.913
%     .911,.727
%     -.966,.788
%     -1.076,1.314
%     -1.524,-.027
%     -1.113,-.030
%     -1.231,.373
%     -.295,-.989
%     .290,-.741
%     .692,-.634
%     1.003,-.553
%     -.893,-1.097
%     -.994,-.856
%     .330,-1.211
%     .211,-.967
%     .466,-.966
%     1.064,-.366
%     .940,-.791];
% 
% orig5_1=[1047,1358];
% z5_1=[870,811
%     1224,949
%     1456,1018
%     1641,1056
%     1444,1370
%     1784,1371
%     1340,1183
%     1448,1137
%     1236,686
%     1319,821
%     1688,1167
%     1072,747
%     1173,825
%     1479,807
%     1595,923
%     458,888
%     400,574
%     126,1376
%     371,1374
%     302,1147
%     866,1947
%     1226,1801
%     1463,1733
%     1649,1692
%     508,2006
%     445,1873
%     1242,2079
%     1173,1931
%     1326,1935
%     1692,1573
%     1608,1828];
% xy5_1=2.54*[z5_1(:,1)-orig5_1(1),orig5_1(2)-z5_1(:,2)]/600;
% xy=xy5_1;
% 
%xy5_2: D=3.4875 cm [angles]=[0 5.45 -1] [0 1 0]
% xy5_2a=2.54*[-.289,.869
%     .297,.660
%     .668,.545
%     .974,.490
%     .659,-.046
%     1.223,-.046
%     .488,.275
%     .671,.347
%     .312,1.090
%     .456,.870
%     1.067,.296
%     .052,.984
%     .205,.853
%     .711,.886
%     .912,.711
%     -.958,.761
%     -1.057,1.255
%     -1.500,-.073
%     -1.108,-.052
%     -1.209,.326
%     -.302,-1.012
%     .307,-.764
%     .716,-.646
%     1.012,-.586
%     -.898,-1.127
%     -.977,-.889
%     .353,-1.242
%     .213,-.983
%     .487,-.990
%     1.076,-.372
%     .955,-.817];
% 
% orig5_2=[1047,1385];
% z5_2=[874,860
%     1225,992
%     1457,1058
%     1642,1095
%     1448,1412
%     1791,1412
%     1342,1224
%     1452,1177
%     1237,733
%     1321,865
%     1693,1209
%     1074,749
%     1171,872
%     1479,854
%     1598,965
%     468,938
%     411,634
%     134,1426
%     378,1421
%     311,1196
%     871,1997
%     1234,1848
%     1472,1778
%     1663,1734
%     506,2062
%     451,1925
%     1256,2131
%     1182,1978
%     1336,1983
%     1704,1614
%     1619,1873];
% xy5_2=2.54*[z5_2(:,1)-orig5_2(1),orig5_2(2)-z5_2(:,2)]/600;
% xy=xy5_2;
%  
% xy5_3:  D=3.4875cm [angles]=[0 5.45 -1.35] [0 1 0]
% orig5_3=[1027,1370];
% z5_3=[857,859
%     1206,988
%     1438,1055
%     1625,1094
%     1429,1408
%     1771,1410
%     1322,1219
%     1431,1174
%     1218,738
%     1303,866
%     1620,1207
%     1056,797
%     1153,872
%     1461,851
%     1578,965
%     453,936
%     0,0
%     118,1420
%     360,1414
%     295,1191
%     848,1995
%     1215,1847
%     1453,1772
%     1642,1736
%     482,2058
%     426,1919
%     1232,2136
%     1158,1979
%     1317,1984
%     1683,1614
%     1602,1876];
% xy5_3=2.54*[z5_3(:,1)-orig5_3(1),orig5_3(2)-z5_3(:,2)]/600;
% xy=xy5_3;
% 
%xy5_4: D=3.4875cm [angles]=[0 5.45 -1.925] [0 1 0]
% orig5_4=[1021,1356];
% z5_4=[853,868
%     1201,995
%     1433,1058
%     1617,1098
%     1425,1409
%     1769,1414
%     1316,1223
%     1426,1177
%     1210,749
%     1294,873
%     1665,1212
%     1051,809
%     1145,879
%     1449,861
%     1572,970
%     455,947
%     408,660
%     110,1427
%     356,1421
%     291,1199
%     844,2008
%     1214,1856
%     1453,1784
%     1644,1740
%     475,2076
%     415,1933
%     1234,2154
%     1161,1992
%     1317,1997
%     1688,1617
%     1606,1890];
% xy5_4=2.54*[z5_4(:,1)-orig5_4(1),orig5_4(2)-z5_4(:,2)]/600;
% xy=xy5_4;
% 
% xy5_5: D=3.44cm [angles]=[.075 5.39 -1.9] [0 1 0]
% orig5_5=[1097,1396];
% z5_5=[930,909
%     1272,1037
%     1504,1100
%     1683,1142
%     1496,1452
%     1836,1457
%     1389,1265
%     1496,1220
%     1286,788
%     1365,915
%     1732,1256
%     1125,848
%     1218,920
%     1519,904
%     1638,1012
%     535,986
%     489,697
%     200,1466
%     441,1460
%     379,1239
%     925,2050
%     1288,1896
%     1527,1825
%     1715,1786
%     558,2119
%     502,1972
%     1310,2197
%     1234,2033
%     1394,2040
%     1753,1661
%     1677,1930];
% xy5_5=2.54*[z5_5(:,1)-orig5_5(1),orig5_5(2)-z5_5(:,2)]/600;
% xy=xy5_5;
% 
% xy5_6:  D=3.44cm [angles]=[.07 4.4 1] [0 1 0]
% orig5_6=[1115,1624];
% z5_6=[973,1016
%     1333,1156
%     1576,1218
%     1767,1256
%     1550,1591
%     1908,1585
%     1444,1400
%     1560,1348
%     0,0
%     1437,1015
%     1811,1379
%     0,0
%     1277,1029
%     1606,998
%     1727,1116
%     574,1103
%     508,774
%     281,1587
%     508,1591
%     435,1367
%     984,2144
%     1330,2014
%     0,0
%     0,0
%     0,0
%     594,2062
%     1348,2272
%     1278,2137
%     1428,2142
%     1805,1790
%     1713,2043];
% xy5_6=2.54*[z5_6(:,1)-orig5_6(1),orig5_6(2)-z5_6(:,2)]/600;
% xy=xy5_6;
% 
% xy5_7: D=3.44cm [angles]=[.07 4.4 1] [0 1 0]
% orig5_7=[1082,1434];
% z5_7=[935,825
%     1297,963
%     1537,1028
%     1733,1064
%     1516,1398
%     1874,1392
%     1411,1208
%     1527,1158
%     0,0
%     1400,823
%     1778,1184
%     0,0
%     1241,836
%     0,0
%     1695,924
%     539,913
%     475,584
%     250,1399
%     474,1401
%     400,1179
%     952,1956
%     1296,1822
%     1534,1758
%     1725,1717
%     620,2000
%     561,1875
%     1316,2082
%     1246,1945
%     1397,1950
%     1771,1597
%     1682,1852];
% xy5_7=2.54*[z5_7(:,1)-orig5_7(1),orig5_7(2)-z5_7(:,2)]/600;
% xy=xy5_7;
% 
% xy5_8: D=3.44cm [angles]=[0 3.41 -3.65] [0 1 0]
% orig5_8=[1086,1439];
% z5_8=[986,1019
%     1321,1131
%     1555,1192
%     1753,1231
%     1559,1547
%     1937,1558
%     1438,1356
%     1553,1311
%     1332,894
%     1415,1013
%     1810,1346
%     1176,955
%     1259,1023
%     1571,998
%     1700,1103
%     618,1104
%     589,846
%     304,1560
%     518,1552
%     468,1344
%     975,2167
%     1351,2013
%     1604,1950
%     1816,1910
%     0,0
%     562,2067
%     0,0
%     1298,2161
%     1468,2179
%     1854,1778
%     1777,2067];
% 
% xy5_8=2.54*[z5_8(:,1)-orig5_8(1),orig5_8(2)-z5_8(:,2)]/600;
% xy=xy5_8;
% 
% xy5_9: D=3.44cm [angles]=[-.05 3.6 -4.5] [0 1 0]
% orig5_9=[1113,1483];
% z5_9=[1010,1093
%     1341,1204
%     1573,1264
%     1771,1306
%     1580,1618
%     1958,1634
%     1464,1425
%     1572,1381
%     1350,978
%     1432,1090
%     1827,1422
%     1195,1033
%     1288,1097
%     1586,1077
%     1712,1181
%     644,1176
%     618,926
%     318,1636
%     536,1625
%     490,1414
%     994,2253
%     1374,2092
%     1631,2026
%     1843,1991
%     0,0
%     573,2157
%     0,0
%     0,0
%     1496,2266
%     1878,1853
%     1813,2154];
% 
% xy5_9=2.54*[z5_9(:,1)-orig5_9(1),orig5_9(2)-z5_9(:,2)]/600;
% xy=xy5_9;
% 
% xy5_10: D=3.44cm [angles]=[-.05 3.775 -5.75] [0 1 0]
% orig5_10=[1093,1450];
% z5_10=[986,1103
%     1311,1212
%     1546,1273
%     1735,1317
%     1556,1625
%     1932,1646
%     1435,1431
%     1544,1388
%     1317,993
%     1401,1101
%     1798,1431
%     1169,1050
%     1258,1107
%     1551,1093
%     1676,1195
%     623,1188
%     600,947
%     282,1651
%     504,1635
%     464,1423
%     962,2295
%     1355,2119
%     1617,2049
%     1829,2013
%     0,0
%     528,2195
%     0,0
%     0,0
%     0,0
%     1861,1867
%     0,0];
% xy5_10=2.54*[z5_10(:,1)-orig5_10(1),orig5_10(2)-z5_10(:,2)]/600;
% xy=xy5_10;
% 
% plot(r(:,1),r(:,2),'.b',xy(:,1),xy(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

%xy7:[angles]=[-1.75,-1.5,-.5] [0 0 1]
% xy7=[-3.266,-.024
%     -2.684,0
%     -3.154,1.094
%     -3.096,-1.110
%     .288,-2.846
%     1.316,-2.920
%     -.714,-2.914
%     -1.256,-2.982
%     1.894,-3.020
%     3.314,.216
%     3.410,-.844
%     3.364,1.300
%     .090,3.086
%     1.130,3.232
%     -.914,3.102
%     1.696,3.390
%     -1.450,3.158
%     .102,2.328
%     .276,-2.120
%     -2.736,-1.512
%     -2.844,1.498];
% plot(r(:,1),r(:,2),'.b',xy7(:,1),xy7(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

%xy9:[angles]=[-3,-2,-1.5] [1 0 0]
% xy9=[-.430,3.616
%     -1.276,3.178
%     -1.756,2.990
%     -2.314,2.760
%     .440,3.206
%     .894,3.026
%     1.430,2.873
%     -.330,2.066
%     .592,2.146
%     -1.292,2.056
%     -2.010,1.094
%     1.314,1.252
%     1.812,1.656
%     -3.166,1.138
%     -2.402,-.316
%     -1.940,-1.750
%     -3.184,-2.020
%     -.094,-2.590
%     1.544,-1.456
%     1.830,-.058
%     -1.848,2.132
%     1.100,2.166
%     -1.780,-2.950
%     1.432,-2.594];
% plot(r(:,1),r(:,2),'.b',xy9(:,1),xy9(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

% xyjc1a1=2.54*[.051,.955
%     -.520,.790
%     -.934,.740
%     -1.289,.728
%     -.567,1.281
%     -.698,1.067
%     -1.542,.188
%     -.996,1.149
%     -.888,.527
%     -.863,.138
%     .589,.781
%     1.023,.702
%     1.390,.664
%     .677,1.259
%     .805,1.025
%     1.627,.088
%     1.115,1.083
%     .975,.495
%     .926,.093
%     -.012,-.681
%     -.519,-.519
%     -.917,-.444
%     -1.245,-.432
%     -.558,-.896
%     -.685,-.709
%     -1.167,-.586
%     .522,-.566
%     .938,-.496
%     1.277,-.460
%     .958,-.810
%     1.179,-.645
%     .557,-.936];
% plot(r(:,1),r(:,2),'.b',xyjc1a1(:,1),xyjc1a1(:,2),'*r')
% axis([-4.4 4.4 -5.7 5.7])

% plot(r(:,1),r(:,2),'.b')
% axis([-4.4 4.4 -5.7 5.7])
%    
% if n_o(1)~=0
%     title('[1 0 0] direction:  c-axis horizontal c=100%, a=100%')
% elseif n_o(2)~=0
%     title('[0 1 0] direction:  c-axis vertical c=100%, a=100%')
% elseif n_o(3)~=0
%     title('[0 0 1] direction c=100%, b=100%, a=100%')
% end
% xlabel('cm')
% ylabel('cm')
% 
% 
% 
% % figure(2)
% % subplot(2,1,1)
% % plot(vbl_dsct,f_Ce,vbl_dsct,f_Ir,vbl_dsct,f_In)
% % title('from reference')
% % legend('Ce','Ir','In')
% % subplot(2,1,2)
% % plot(data(:,9),f1,'.b',data(:,9),f2,'.r',data(:,9),f3,'k.')
% % title('points for our particular poles')
% % legend('Ce','Ir','In')
% 
% 
% display('1     2      3       4      5       6      7       8             9                10           11        12        13')
% display('m1    m2     m3     [nx  ,  ny  ,  nz]   theta   lambda   sin(theta)/lambda   magnitude_n   Intensity  x-coord.  y-coord.');
