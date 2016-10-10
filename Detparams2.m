clc; 
clear all; 
close all;

pi=3.14159265358979;

%% POLE LOCATIONS 

%%SAMPLE 5
%%5_1
orig5_1=[1047,1358];
z5_1=[870,811
    1224,949
    1456,1018
    1641,1056
    1444,1370
    1784,1371
    1340,1183
    1448,1137
    1236,686
    1319,821
    1688,1167
    1072,747
    1173,825
    1479,807
    1595,923
    458,888
    400,574
    126,1376
    371,1374
    302,1147
    866,1947
    1226,1801
    1463,1733
    1649,1692
    508,2006
    445,1873
    1242,2079
    1173,1931
    1326,1935
    1692,1573
    1608,1828];
xy5_1=2.54*[z5_1(:,1)-orig5_1(1),orig5_1(2)-z5_1(:,2)]/600;
xy=xy5_1;
zxy=z5_1;
%%%CONVERGED
n_owe=-1*[0 1 0];
alph=(pi/180)*.016344;
bet=(pi/180)*5.6324;
gam=(pi/180)*-.47451;
D=3.4934*10^(-2);
g=12;

%5_10

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
% zxy=z5_10;
% xy=xy5_10;
% % % CONVERGED
% n_owe=-1*[0 1 0];
% alph=(pi/180)*-.18563;
% bet=(pi/180)*3.7577;
% gam=(pi/180)*-5.7434;
% D=3.4202*10^(-2);
% g=12;


%% HKL VALUES 

%Sample 5
hkl=[0,2,1
        1,5,2
        1,3,1
        3,7,2
        1,3,0
        1,2,0
        3,11,2
        3,9,2
        1,5,3
        1,4,2
        5,11,2
        1,9,5
        1,6,3
        2,6,3
        2,5,2
        -1,5,2
        -1,5,3
        -1,3,0
        -1,4,0
        -3,11,2
        0,2,-1
        1,5,-2
        1,3,-1
        3,7,-2
        -1,6,-3
        -1,5,-2
        1,5,-3
        1,6,-3
        1,4,-2
        5,11,-2
        2,5,-2];
 n=[1,2,3,4,5,6];
 
eps=2;
overall=0;
converge=0;
conv=0;
oall=0;

%% TO CHECK HOW/WHEN THE PARAMETERS CONVERGE

exeunt=0;
exita=0;
exitb=0;
exitg=0;
exitD=0;
where0=0;
where1=0;
where2=0;
where3=0;

%CHOOSE INCIDENT BEAM DIRECTION

% When analyzing the pressure effects, graphing with the c-axis into the
% plane of the image, the "a" dimension runs vertically and the "b"
% dimension runs horizontally.

% % % dink=1;
% % % while dink==1
% % %     n_owe=-1*input('input the incident beam direction: only either [1,0,0] , [0,1,0] , [0,0,1] allowed:  ');
% % %     
% % % %     [1 0 0] has c-axis running horizontally, [0 1 0] c-axis runs vertically
% % %     
% % %     if sum(n_owe)~=-1 | (n_owe*n_owe.')~=1
% % %         display('try again');
% % %     elseif n_owe(1)==-1 & n_owe(2)==0 & n_owe(3)==0
% % %         dink=2;
% % %     elseif n_owe(1)==0 & n_owe(2)==-1 & n_owe(3)==0
% % %         dink=2;
% % %     elseif n_owe(1)==0 & n_owe(2)==0 & n_owe(3)==-1
% % %         dink=2;
% % %     else
% % %         display('try again');
% % %     end
% % % end

%INPUT ANGLES

% % % alph=(pi/180)*input('input angle alpha in degrees, see code for angle definition:  ');
% % % bet=(pi/180)*input('input angle beta in degrees:  ');
% % % gam=(pi/180)*input('input angle gamma in degrees:  ');
%The angle alpha is the rotation about the incident beam direction such
%that, as viewed from that direction.  Beta is the rotation about the crystal axis aligned with the "y" axis on the image when alpha =0.  
%Gamma is the rotation about the crystal axis aligned with the "x" axis on the image when alpha =0.  Of course the order of
%operation of the rotation matrices affects the resultant crystal
%orientation.  They are set to operate in a particular order:  1st=Gamma,
%2nd=Beta, 3rd=Alpha.  However, one nonzero angle at a time is best and the
%rotation becomes the identity matrix when an angle is zero.

% max n_o||x,y is 26; n_o || z is 42
% % % g=input('type the abs(maximum) hkl value you wish to be evaluated;optimum=[0 0 1]-->30,[0 1 0] | [1 0 0]-->20:  ');

%% INITIALIZE STUFF 1

alpha=alph;
betb=bet;
gamg=gam;
DD=D;
D1=0;
D2=0;

%% INITIAL VALUE WHEN TO EXIT

lilang=pi/180*10^(-4);
lildist=10^(-6);
lilbig=10^(-7);
delta=0;
deltb=0;
deltg=0;
deltD=0;
stepDD=.05*DD;

%% INITIALIZE PARAMETERS THAT INCREASE BY OVERALL OR CONVERGE

maxo=400;
params=zeros(maxo,9);
bigdiff=zeros(maxo-1,1);
mbigdiff=zeros(maxo-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while eps>1
    overall=overall+1
    data=press115(alph,bet,gam,D,n_owe,g);

    if overall==1
        figure(1)
        plot(data(:,12),data(:,13),'.b',xy(:,1),xy(:,2),'*r')
        axis([-4.4 4.4 -5.7 5.7])
        title('[1 0 0] c-axis horizontal:  first guess')
        xlabel('cm')
        ylabel('cm')
        pause(5)
    end


    %% MATCHING EXPERIMENTAL AND SIMULATED POINTS BY HKL VALUE
    %% DEFINING VECT 
    cnt=zeros(1,size(hkl,1));
    zers=0;
    midpoint=0;
    vect=[];
    vect=zeros(size(hkl,1),19);
    clump=0;

    for j=1:size(hkl,1); % 1 -- scan thru experimental points, hkl value
        for k=1:length(n); % 2 -- scan thru orders of constructive interference
            for i=1:size(data,1); % 3 -- scan thru data file (containing ALL simulated points)

                if n(k)*hkl(j,1)==data(i,1) & n(k)*hkl(j,2)==data(i,2) & n(k)*hkl(j,3)==data(i,3)
                    clump=clump+1;
                    if zxy(j,1)==0 && zxy(j,2)==0  %% ELIMINATING THE HKL POINTS WITH NO EXPERIMENTAL VALUE
                        zers=zers+1;
                        cnt(j)=1;
                        vect(j,:)=[zxy(j,1),zxy(j,2),n(k),hkl(j,:),data(i,:)]; % Sets vect equal to zero
                        break % Break out of scanning thru data files (3) AND orders, n (2)
                    else
                        cnt(j)=1;
                        vect(j,:)=[xy(j,1)-data(i,12),xy(j,2)-data(i,13),n(k),hkl(j,:),data(i,:)];  % Contains 19 columns
                        break % Break out of scanning thru data files (3) and orders, n (2)
                    end
                end
                if midpoint==0;
                    if data(i,1)~=0 && data(i,2)==0 && data(i,3)==0;
                        midpoint=1;
                        nill(1,:)=data(i,:);
                    end
                end                  

            end % 3

            if cnt(j)==1;
                break
            end

        end % 2

    end % 1

    tic1=1;

    %% INITIALIZE SIZE OF GROUP (FIRST INDEX # DEPENDS ON WHICH POLES ARE MARKED)

    %%Sample 5
    midy=zeros(4,19);
    diff=zeros(1,size(vect,1));
    midiff=zeros(1,size(midy,1));

    for j1=1:size(vect,1);
        %% Separate y-axis hkl values of zero  (ALWAYS y-axis: hkl is ORIENTATION
        %% DEPENDENT)

        if n_owe(1)~=0
            %%(k=0)

            if vect(j1,8)==0
                midy(tic1,:)=vect(j1,:);
                midiff(tic1)=sqrt(midy(tic1,1)^2+midy(tic1,2)^2);
                tic1=tic1+1;
            end

        elseif n_owe(2)~=0
            %%(l=0)

            if vect(j1,9)==0
                midy(tic1,:)=vect(j1,:);
                midiff(tic1)=sqrt(midy(tic1,1)^2+midy(tic1,2)^2);
                tic1=tic1+1;
            end

        elseif n_owe(3)~=0
            %%(h=0)
            
            if vect(j1,7)==0
                midy(tic1,:)=vect(j1,:);
                midiff(tic1)=sqrt(midy(tic1,1)^2+midy(tic1,2)^2);
                tic1=tic1+1;
            end
        end

        diff(j1)=sqrt(vect(j1,1)^2+vect(j1,2)^2);

    end

    %% CALCULATE THE CONVERGENCE TEST FUNCTION
    
    bigdiff(overall)=sum(diff);
    mbigdiff(overall)=sum(midiff);
    
    %% TEST FOR CONVERGENCE
            
    if overall-converge~=1
        
        if exitg==0
            if mod(overall-converge,3)==1
                if abs((bigdiff(overall-3)-bigdiff(overall)))/bigdiff(overall-3)<=lilbig
                    exeunt=exeunt+1;
                    if exeunt==1
                        exita=1
                        alph*(180/pi)
                    elseif exeunt==2
                        exitb=1
                        bet*(180/pi)
                    elseif exeunt==3
                        exitg=1
                        gam*(180/pi)
                        
                        oall=overall;
                        
                    end
                end
            end
        elseif exitg~=0 && exitD==0 && overall-oall>3
            if abs(D1-D2)<lildist && converge==0
                conv=conv+1
                D=.5*(D1+D2);
                D*100
                exeunt=0;
                exita=0;
                exitb=0;
                exitg=0;
                exitD=0;
                oall=0;
                stepDD=stepDD/(conv+1);
                converge=overall-1;
            elseif abs(D1-D2)<lildist && converge~=0
                D=.5*(D1+D2);
                D*100
                if abs(params(converge,1)-params(overall,1))<=(180/pi)*lilang && abs(params(converge+1,2)-params(overall,2))<=(180/pi)*lilang && abs(params(converge+1,3)-params(overall,3))<=(180/pi)*lilang && abs(params(converge+1,4)-100*D)<=100*lildist
                    display('CONVERGED')
                    exitD=1;
                    eps=.5;
                    conv=conv+1
                else
                    conv=conv+1
                    stepDD=stepDD/(conv+1);
                    exeunt=0;
                    exita=0;
                    exitb=0;
                    exitg=0;
                    exitD=0;
                    oall=0;
                    converge=overall-1;
                end
            
            end
        end
    end
       
    %% CHANGE PARAMETERS
    
    if exitg==0
        if mod(overall-converge,3)==1
            
            % x+dx    

            if exita==0

                delta=lilang;
                alph2=alph;
                alph=alph+delta;

            elseif exita~=0 && exitb==0

                deltb=lilang;
                bet2=bet;
                bet=bet+deltb;

            elseif exitb~=0 && exitg==0

                deltg=lilang;
                gam2=gam;
                gam=gam+deltg;
            end

        elseif mod(overall-converge,3)==2

            % x-dx

            if exita==0

                delta=-lilang;
                alph=alph2+delta;

            elseif exita~=0 && exitb==0

                deltb=-lilang;
                bet=bet2+deltb;

            elseif exitb~=0 && exitg==0

                deltg=-lilang;
                gam=gam2+deltg;
            end

        elseif mod(overall-converge,3)==0

            % x+correction

            if exita==0
                fprime=(bigdiff(overall-1)-bigdiff(overall-2));
                f2prime=(bigdiff(overall-1)-2*bigdiff(overall-2)+bigdiff(overall))/lilang;
                delta=-fprime/f2prime;
                alph=alph2+delta;

            elseif exita~=0 && exitb==0
                fprime=(bigdiff(overall-1)-bigdiff(overall-2));
                f2prime=(bigdiff(overall-1)-2*bigdiff(overall-2)+bigdiff(overall))/lilang;
                deltb=-fprime/f2prime;
                bet=bet2+deltb;

            elseif exitb~=0 && exitg==0
                fprime=(bigdiff(overall-1)-bigdiff(overall-2));
                f2prime=(bigdiff(overall-1)-2*bigdiff(overall-2)+bigdiff(overall))/lilang;
                deltg=-fprime/f2prime;
                gam=gam2+deltg;
            end

        end
        
    elseif exitg~=0 && exitD==0
              
        if overall-oall<=1 %Establish the left boundary - D1
            if mod(overall-oall,2)==0
                deltD=lildist;
                D1=D;
                D=D+deltD;
            elseif mod(overall-oall,2)==1 %Use slope at left boundary to pick direction to alter D
                fprimeD=(mbigdiff(overall)-mbigdiff(overall-1));
                deltD=-sign(fprimeD)*stepDD; % Initial stepD defined above - stepD=.05*DD -- % Step towards the downslope
                D=D1+deltD;
                D2=D;
                stepD=stepDD;%Establish evaluation boundary - D2 (not necessarily to the right of D1)
                
            end
                
        else
            if mod(overall-oall,2)==0 %Find slope at D2
                deltD=-sign(fprimeD)*lildist;
                D=D+deltD;
            elseif mod(overall-oall,2)==1
                if (mbigdiff(overall)-mbigdiff(overall-1))/lildist<-lilbig % Slope at D2 is the same as at D1, initially -> shift in that direction
                    D1=D2;
                    deltD=-sign(fprimeD)*stepD;
                    D2=D1+deltD;
                    D=D2;

                    where1=1;
                    
                elseif (mbigdiff(overall)-mbigdiff(overall-1))/lildist>lilbig % Slope at D2 is opposite D1 -> cut step in half
                    stepD=stepD/2;
                    deltD=-sign(fprimeD)*stepD;
                    D2=D1+deltD;
                    D=D2;
                    
                    where2=1;
                    
                else
                    display('converged within error?')
                    
                    where3=1;
                                        
                    figure(2)
                    plot(data(:,12),data(:,13),'.b',xy(:,1),xy(:,2),'*r')
                    axis([-4.4 4.4 -5.7 5.7])
                    title('converged within error?')
                    pause
                end
            end
        end
    end

    params(1,:)=[alpha*180/pi betb*180/pi gamg*180/pi DD*100 0 0 0 0 0];
    params(overall+1,:)=[alph*180/pi bet*180/pi gam*180/pi D*100 D1*100 D2*100 where1 where2 where3];
    where1=0;
    where2=0;
    where3=0;
    
    if overall==maxo
        figure(2)
        plot(data(:,12),data(:,13),'.b',xy(:,1),xy(:,2),'*r')
        axis([-4.4 4.4 -5.7 5.7])
        if n_owe(1)~=0
            title('[1 0 0] c-axis horizontal:  maxed')
        elseif n_owe(2)~=0
            title('[0 1 0] c-axis vertical:  maxed')
        elseif n_owe(3)~=0
            title('[0 0 1] a-axis vertical, b-axis horizontal:  maxed')
        end
        xlabel('cm')
        ylabel('cm')
        return
    end
end
figure(2)
plot(data(:,12),data(:,13),'.b',xy(:,1),xy(:,2),'*r')
axis([-4.4 4.4 -5.7 5.7])
if n_owe(1)~=0
    title('[1 0 0] c-axis horizontal:  converged')
elseif n_owe(2)~=0
    title('[0 1 0] c-axis vertical:  converged')
elseif n_owe(3)~=0
    title('[0 0 1] a-axis vertical, b-axis horizontal:  converged')
end
xlabel('cm')
ylabel('cm')