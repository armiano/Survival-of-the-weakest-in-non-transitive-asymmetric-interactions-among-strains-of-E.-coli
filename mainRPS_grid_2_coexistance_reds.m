
%%  
%set 3 grids, one for each strain. It will be binary, so either the strain
%is present (1) or not (0)
close all
clear all
N = 150;
gridR = zeros(N,N);
gridP = zeros(N,N);
gridS = zeros(N,N);
PLOT =0;
SAVE =0;
VIDEO = 0;
Ntrial = 100;

nx = N ; dx = 1;
ny = N ; dy = 1;

xgrid = 1:dx:nx;
ygrid = 1:dy:ny;
[X, Y] = meshgrid(xgrid);

radius = 25; %range to try from 6 to 23 %for the version without (1-n1) is from to 16
radius2 = 50;
radius3 = 75;
distance = 0;
positiony = ny/2;
positionx1 = nx/2 + distance;
positionx2 = nx/2 - distance;
for i = 1:nx
    for j = 1:ny
          
        if ((X(i,j)-positionx1).^2+(Y(i,j)-positiony).^2) < radius^2 
            gridR(i,j) = 1;
            
        end
        
        if ((X(i,j)-positionx1).^2+(Y(i,j)-positiony).^2) < radius2^2 && ((X(i,j)-positionx1).^2+(Y(i,j)-positiony).^2) > radius^2
            gridS(i,j) = 1;
            
        end
        
         if ((X(i,j)-positionx1).^2+(Y(i,j)-positiony).^2) < radius3^2 && ((X(i,j)-positionx1).^2+(Y(i,j)-positiony).^2) > radius2^2
            gridP(i,j) = 1;
            
        end
    end
end

%figure;imagesc(gridR); drawnow

%set 0 BC
gridR(1,:)=0;
gridR(N,:)=0;
gridR(:,1)=0;
gridR(:,N)=0;

gridP(1,:)=0;
gridP(N,:)=0;
gridP(:,1)=0;
gridP(:,N)=0;

gridS(1,:)=0;
gridS(N,:)=0;
gridS(:,1)=0;
gridS(:,N)=0;

%plot them all together
figure
rgbImage = cat(3, gridR, gridS,gridP);  %need to get it as a 3D matrix to get them to overlap
imagesc(rgbImage);
title('Together')



gridRinit = gridR;
gridPinit = gridP;
gridSinit = gridS;



%%
%we keep pR and pP which correspond to the strongest and the weakest 
pR =  linspace(0.05,0.1,8);%probability of strain R to be killed by P
pP =  linspace(0.05,0.1,8); %probability of strain P to be killed by S
pS =  linspace(0.05,0.5,8);%probability of strain S to be killed by R

writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj);

arraypP = 0.1:0.05:0.5;
%arraypP=0.4;
%we will sweep across all values between the range

steadystate = zeros(length(arraypP),Ntrial,3);

%% Run Loop
%need to go through the entire maps and check each pixel for each of the 3
%matrix


for arr = 1:length(arraypP)

pP =  linspace(0.05,arraypP(arr),8);


for trial = 1:Ntrial
    
gridR=gridRinit; %reset initial condition array
gridP=gridPinit;
gridS=gridSinit;
    
sumGrid = gridR + gridP +gridS; %check if there is an empty spot in the overall grid

  disp(['pP',num2str(arraypP(arr)), '_trial', num2str(trial)])
for t = 1:10000
    
        R = gridR;
        P = gridP;
        S = gridS;
        sumGrid = R + P +S;
        gridRTotFrac(t) = sum(sum(R))/(N^2);
        gridPTotFrac(t) = sum(sum(P))/(N^2);
        gridSTotFrac(t) = sum(sum(S))/(N^2);
        
    for a = 2:N-1
        for b = 2:N-1
            

            if sumGrid(a,b) == 0  %if the pixel is empty, it is populated according to a set of weighted probabilities proportional to the occupancy of each strain at the eight neighboring locations
               sumR = R(a-1,b-1)+R(a-1,b)+R(a,b-1)+R(a+1,b+1)+R(a,b+1)+R(a+1,b)+R(a+1,b-1)+R(a-1,b+1); % sum of 8 surrounding pixels
               sumP = P(a-1,b-1)+P(a-1,b)+P(a,b-1)+P(a+1,b+1)+P(a,b+1)+P(a+1,b)+P(a+1,b-1)+P(a-1,b+1);
               sumS = S(a-1,b-1)+S(a-1,b)+S(a,b-1)+S(a+1,b+1)+S(a,b+1)+S(a+1,b)+S(a+1,b-1)+S(a-1,b+1);
               
    %            
               pRFill =  sumR/8;%probability of being filled with R
               pPFill =  sumP/8; %probability of being filled with P
               pSFill =  sumS/8; %probability of being filled with S

              % gives back a number according to the probablities 1=R, 2=P, 3=S 
               if pRFill==0 && pPFill == 0 && pSFill ==0
                   y = 0;
                  % disp(y)
               else
                    y = randsample([1 2 3],1,true,[pRFill pPFill pSFill]);
               end

               if y==1
                   gridR(a,b) = 1;
               end

               if y == 2
                   gridP(a,b) = 1;
               end

               if y == 3
                   gridS(a,b) = 1;
               end
               
             
            else %if the spot is not empty, one of the matrix must be 1
                
                if R(a,b) ==1
                  sumP = P(a-1,b-1)+P(a-1,b)+P(a,b-1)+P(a+1,b+1)+P(a,b+1)+P(a+1,b)+P(a+1,b-1)+P(a-1,b+1);
                 if sumP == 0
                     Prob = 0.05; %basal probability
                 else
                     if sumP > 4
                     Prob = pR(8);
                     else
                     Prob = pR(4);
                     end
                 end
                  x = rand; %rand variable between 0 and 1
                  if x < Prob
                      gridR(a,b) =0; %killed
                  end
               
                  
                  
                elseif P(a,b) ==1
                   sumS = S(a-1,b-1)+S(a-1,b)+S(a,b-1)+S(a+1,b+1)+S(a,b+1)+S(a+1,b)+S(a+1,b-1)+S(a-1,b+1);
                 if sumS == 0
                     Prob = 0.05; %basal probability
                 else
                     if sumS > 4
                     Prob = pP(8);
                     else
                     Prob = pP(4);
                     end
                 end
                  x = rand; %rand variable between 0 and 1
                  if x < Prob
                      gridP(a,b) =0; %killed
                  end
                    
                  
                  
                  
                elseif S(a,b) ==1 
                   sumR = R(a-1,b-1)+R(a-1,b)+R(a,b-1)+R(a+1,b+1)+R(a,b+1)+R(a+1,b)+R(a+1,b-1)+R(a-1,b+1); % sum of 8 surrounding pixels
                 if sumR == 0
                     Prob = 0.05; %basal probability
                 else
                     if sumR > 4
                     Prob = pS(8);
                     else
                     Prob = pS(4);
                     end
                 end
                  x = rand; %rand variable between 0 and 1
                  if x < Prob
                      gridS(a,b) =0; %killed
                  end
                    
                end
             

            

            end
        end
    end
    
    
   
    
end

 if SAVE ==1
    figure
    plot(gridRTotFrac,'r')
    hold on
    plot(gridPTotFrac,'b')
    hold on 
    plot(gridSTotFrac,'g')
    hold off
    path=['/Users/ariannamiano/Desktop/PhD/ModelRPSPlate/Coexistance_sweep_circle/pR_0.1_pS_0.5_pP_from_0.4_to_0.5/pP' num2str(arraypP(arr)) '_trial' num2str(trial) '.svg'] ;
    saveas(gcf,path);
    close all
 end   
 
 if PLOT==1
     figure
    plot(gridRTotFrac,'r')
    hold on
    plot(gridPTotFrac,'b')
    hold on 
    plot(gridSTotFrac,'g')
    hold off
 end
    steadystate(arr,trial,1) = gridRTotFrac(end);
    steadystate(arr,trial,2) = gridPTotFrac(end);
    steadystate(arr,trial,3) = gridSTotFrac(end);
    %this matrix will store the last value of each species for each
    %parameteer pP and for each trial
  
end
end

figure
plot(steadystate(2,:,2),'b')
hold on
plot(steadystate(2,:,1),'r')
hold on
plot(steadystate(2,:,3),'g')
title(num2str(arraypP(2)))

figure
plot(steadystate(1,:,2),'b')
hold on
plot(steadystate(1,:,1),'r')
hold on
plot(steadystate(1,:,3),'g')
title(num2str(arraypP(1)))

figure
plot(steadystate(7,:,2),'b')
hold on
plot(steadystate(7,:,1),'r')
hold on
plot(steadystate(7,:,3),'g')
title(num2str(arraypP(7)))


