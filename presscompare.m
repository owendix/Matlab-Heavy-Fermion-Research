clc;clear all;close all;

pi=3.14159265358979;
Dnot=3.45*10^(-2);

%ROTATES EXPERIMENTAL POINTS FOR VISUAL COMPARISON OF LATTICE CONSTANT
%SHIFTS

for overall=1:2
    if overall==1
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
        n_o=[0 -1 0];
        
        % ANGLES ARE NEGATIVE OF THEIR CONVERGED VALUE FROM Det_params.m
        alph=-(pi/180)*.014596;
        bet=-(pi/180)*5.6324;
        gam=(pi/180)*.47521;
        D=3.4934*10^(-2);
        
        alph=-(pi/180)*.016344;
        bet=-(pi/180)*5.6324;
        gam=(pi/180)*.47451;
        D=3.4934*10^(-2);
        
        
    elseif overall==2
        
        orig5_10=[1093,1450];
        z5_10=[986,1103
            1311,1212
            1546,1273
            1735,1317
            1556,1625
            1932,1646
            1435,1431
            1544,1388
            1317,993
            1401,1101
            1798,1431
            1169,1050
            1258,1107
            1551,1093
            1676,1195
            623,1188
            600,947
            282,1651
            504,1635
            464,1423
            962,2295
            1355,2119
            1617,2049
            1829,2013
            0,0
            528,2195
            0,0
            0,0
            0,0
            1861,1867
            0,0];
        xy5_10=2.54*[z5_10(:,1)-orig5_10(1),orig5_10(2)-z5_10(:,2)]/600;
        xy=xy5_10;
        %CONVERGED
        n_o=[0 -1 0];
        alph=(pi/180)*.18820;
        bet=-(pi/180)*3.7473;
        gam=(pi/180)*5.7398;
        D=3.4232*10^(-2);
    end
    %ROTATION MATRICES DEFINED


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

    %TRANSFORM TO 3D, ROTATE, TRANSFORM TO 2D

    if n_o(1)~=0
        ka=0;
        for k=1:size(xy,1)
            if xy(k,1)~=0 || xy(k,2)~=0
                ka=ka+1;
                n(:,ka)=(1/sqrt((100*D)^2+xy(k,1)^2+xy(k,2)^2))*[100*D;xy(k,2);xy(k,1)]; %TRANSFORM TO 3D
                enn(:,ka)=n(:,ka)-n_o.';
                enrot(:,ka)=Rx*Ry*Ro*enn(:,ka); %ROTATE
                nrot(:,ka)=enrot(:,ka)+n_o.';
                xyrot(ka,:)=100*Dnot*[nrot(3,ka)./abs(nrot(1,ka)),nrot(2,ka)./abs(nrot(1,ka))]; %TRANSFORM TO 2D
            end
        end
    elseif n_o(2)~=0
        ka=0;
        for k=1:size(xy,1)
            if xy(k,1)~=0 || xy(k,2)~=0
                ka=ka+1;
                n(:,ka)=(1/sqrt((100*D)^2+xy(k,1)^2+xy(k,2)^2))*[xy(k,1);100*D;xy(k,2)]; %TRANSFORM TO 3D
                enn(:,ka)=n(:,ka)-n_o.';
                enrot(:,ka)=Rx*Ry*Ro*enn(:,ka); %ROTATE
                nrot(:,ka)=enrot(:,ka)+n_o.';
                xyrot(ka,:)=100*Dnot*[nrot(1,ka)./abs(nrot(2,ka)),nrot(3,ka)./abs(nrot(2,ka))]; %TRANSFORM TO 2D
            end
        end
    elseif n_o(3)~=0
        ka=0;
        for k=1:size(ord_data,1)
            if xy(k,1)~=0 || xy(k,2)~=0
                ka=ka+1;
                n(:,ka)=(1/sqrt((100*D)^2+xy(k,1)^2+xy(k,2)^2))*[xy(k,2);xy(k,1);100*D]; %TRANSFORM TO 3D
                enn(:,ka)=n(:,ka)-n_o.';
                enrot(:,ka)=Rx*Ry*Ro*enn(:,ka); %ROTATE
                nrot(:,ka)=enrot(:,ka)+n_o.';
                xyrot(ka,:)=100*Dnot*[nrot(2,ka)./abs(nrot(3,ka)),nrot(1,ka)./abs(nrot(3,ka))]; %TRANSFORM TO 2D
            end
        end
    end
    
    if overall==1
        xyrot1=xyrot;
        nrot1=nrot;
        enrot1=enrot;
        enn1=enn;
        n1=n;
    elseif overall==2
        xyrot2=xyrot;
        nrot2=nrot;
        enrot2=enrot;
        enn2=enn;
        n2=n;
        alph=0;
        bet=0;
        gam=0;
        D=Dnot;
        n_owe=[0 -1 0];
        g=12;
        data=press115(alph,bet,gam,D,n_owe,g);
    end        
end
plot(data(:,12),data(:,13),'.b',xyrot1(:,1),xyrot1(:,2),'*g',xyrot2(:,1),xyrot2(:,2),'*r')
axis([-4.4 4.4 -5.7 5.7])
if n_o(1)~=0
    title('[1 0 0] direction:  c-axis horizontal')
elseif n_o(2)~=0
    title('[0 1 0] direction:  c-axis vertical')
elseif n_o(3)~=0
    title('[0 0 1] direction')
end
xlabel('cm')
ylabel('cm')
legend('simulated','P~0','P~5.3kbar')
