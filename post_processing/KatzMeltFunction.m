%crude script to construct melting function of Katz et al. 2003
%pressure is in GPa , Temperature in C, includes parameterization of
%hydrous melting through xH2O

%produces a matrix of melt fraction values at specified pressure and
%temperature, water content is a parameter

clear all;
close all;

alg='trust-region-dogleg';%{'levenberg-marquardt',.01};%'trust-region-dogleg'

options = optimset('algorithm',alg,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1000);

%weight pct water
xH2O = 0.1;

numP=42;
numT=22;

F=zeros(numP,numT);
Fcpx=zeros(numP,numT);
Press=zeros(numP,numT);
Temp=zeros(numP,numT);

dP = 0.05;

MinP = 0.2; %GPa
MaxP = 10;
MinT = 850;
MaxT = 1700; %C

Plist = linspace(MinP,MaxP,numP-2); %pressure steps in GPa
Tlist = linspace(MinT,MaxT,numT-2); %temperature steps in Celsius

dp = Plist(2)-Plist(1);
dT = Tlist(2)-Tlist(1);
Plist = [Plist(1)-dp Plist Plist(end)+dp];
Tlist = [Tlist(1)-dT Tlist Tlist(end)+dT];


flast = 0;
for ii=1:numP
    for jj=1:numT
        
        P=Plist(ii);
        T=Tlist(jj);
        
        M = .15;
        r0 = 0.5;
        r1 = .08;
        
        
        %solidus
        %A1 = 1085.7;
        %A2 = 132.9;
        %A3 = -5.1;
        %k = 43;
        %gamma = .75;
        %DH2O = 0.01;
        %XH2Obulk=0.00;
        
        Fcpxout = M/(r0 + r1*P);
        ff=@(x) f1(x,P,T,xH2O); %cpx function
        
        [MinCpx,X,exitflag] =fsolve(ff,1e-5,options);
        
        if(real(MinCpx)>=0)
            if(MinCpx<=Fcpxout)
                F(ii,jj)=real(MinCpx);
                Fcpx(ii,jj)=real(MinCpx);
            else
                ff=@(x) f2(x,P,T,xH2O); %opx function
                MinOpx=fsolve(ff,flast,options);
                if(imag(MinOpx)>1e-10)
                    disp('Imag value, solver fails to meet tolerance. Likely near cpx-out kink in productivity.')
                end
                if(MinOpx<1)
                    F(ii,jj)=real(MinOpx); %this is the melt fraction
                else
                    F(ii,jj)=1;
                end
            end
        else
            F(ii,jj)=0;
        end
        % increment P slightly
        P = P + dP;
        ff=@(x) f1(x,P,T,xH2O); %cpx function
        [MinCpx,X,exitflag] =fsolve(ff,0,options);
        if(real(MinCpx)>=0)
            if(MinCpx<=Fcpxout)
                fp=real(MinCpx);
%                 Fcpx(ii,jj)=real(MinCpx);
            else
                ff=@(x) f2(x,P,T,xH2O); %opx function
                MinOpx=fsolve(ff,flast,options);
                if(imag(MinOpx)>1e-10)
                    disp('Imag value, solver fails to meet tolerance. Likely near cpx-out kink in productivity.')
                end
                if(MinOpx<1)
                    fp=real(MinOpx); %this is the melt fraction
                else
                    fp=1;
                end
            end
        else
            fp=0;
        end
        dFdp(ii,jj) = (fp-F(ii,jj))/dP;
        
        Press(ii,jj)=P;
        Temp(ii,jj)=T;
        flast = F(ii,jj);
        %Tsolidus(ii,jj) = A1 + A2*P + A3*P^2 - k*(XH2Obulk/(DH2O ))^gamma;
            
        end
    end
    figure;
    surf(Press,Temp,F); xlabel('Pressure (GPa)'); ylabel('Temperature (C)'); zlabel('Melt fraction F')
    
    %TotMelt(:,1)=Press(F(:)>0&F(:)<1);
    %TotMelt(:,2)=Temp(F(:)>0&F(:)<1);
    %TotMelt(:,3)=F(F(:)>0&F(:)<1);
    
    melt_table.P = Press(2:end-1,2:end-1);
    melt_table.T = Temp(2:end-1,2:end-1);
    melt_table.F = F(2:end-1,2:end-1);
    dFdp1 = zeros(numP-2,numT-2);
    dFdT1 = zeros(numP-2,numT-2);
    for ii=2:numP-1
        for jj=2:numT-1
            dFdp1(ii-1,jj-1) = (F(ii+1,jj)-F(ii-1,jj))/(Plist(ii+1)-Plist(ii-1));
            dFdT1(ii-1,jj-1) = (F(ii,jj+1)-F(ii,jj-1))/(Tlist(jj+1)-Tlist(jj-1));
        end
    end
    melt_table.dFdP = dFdp1;
    melt_table.dFdT = dFdT1;
    save(['melt_table_' num2str(xH2O) '.mat'],'melt_table');
    
    
