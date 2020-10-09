%%  
%set 3 grids, one for each strain. It will be binary, so either the strain
%is present (1) or not (0)
close all
clear all
N = 150;
gridR = zeros(N,N);
gridP = zeros(N,N);
gridS = zeros(N,N);
PLOT =1;
SAVE =0;
VIDEO = 0;

gridR(1:75,1:75) = 1;
gridP(76:150,1:75)= 1;
gridS(37:112,76:150)= 1;



%plot them all together
figure
rgbImage = cat(3, gridR, gridS,gridP);  %need to get it as a 3D matrix to get them to overlap
imagesc(rgbImage);
title('Together')

%%
%R is red
%P is blue
%S is green

%define probability of killing, they are the max value when the cell is
%fully surrounded
%lower value correspond to the probability of dying if no killers are
%around
%top value is the highest probability if all cells are around
pR =  linspace(0.05,0.1,8);%probability of strain R to be killed by P
pP =  linspace(0.05,0.28,8); %probability of strain P to be killed by S
pS =  linspace(0.05,0.43,8);%probability of strain S to be killed by R
 
 %we will define toxicity as a parameter whose values is defined by number
 %of neighbours around. The max number 8 means certain death (probability
 %of 1). Therefore toxicity will be defined individually for each cell type
 %as an equally spaced array from 0.1/p (basal death) to 1/p. 
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj);

%% Run Loop
%need to go through the entire maps and check each pixel for each of the 3
%matrix

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
% % 
%  gridR(1,:)=gridR(N,:);
%  gridR(:,N)=gridR(:,1);
%  gridP(1,:)=gridP(N,:);
%  gridP(:,N)=gridP(:,1);
%  gridS(1,:)=gridS(N,:);
%  gridS(:,N)=gridS(:,1);

update = zeros(N,N);
n = (N^2);  % choose the amount of cells to be updated (a third of the total)
update(randperm(numel(update), n)) = 1;
ind =0;
sumGrid = gridR + gridP +gridS; %check if there is an empty spot in the overall grid
for t = 1:2000
    
        R = gridR;
        P = gridP;
        S = gridS;
        sumGrid = R + P +S;
        gridRTotFrac(t) = sum(sum(R))/(N^2);
        gridPTotFrac(t) = sum(sum(P))/(N^2);
        gridSTotFrac(t) = sum(sum(S))/(N^2);
        
    for a = 2:N-1
        for b = 2:N-1
            
            if update(a,b) == 1 %only update in the positions randomly selected

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
    
    if mod(t,5) ==1 && PLOT == 1
    rgbImage = cat(3, gridR,gridS,gridP);  %need to get it as a 3D matrix to get them to overlap
    imagesc(rgbImage);
    title(['time:',num2str(t),'Red=',num2str(round(gridRTotFrac(t),2)),'%, Blue=',num2str(round(gridPTotFrac(t),2)),'%, Green=',num2str(round(gridSTotFrac(t),2)),'%'])
     drawnow
    pause(0.2)
    
    elseif mod(t,10) ==1 && SAVE == 1
    rgbImage = cat(3, gridR,gridS,gridP);
    imagesc(rgbImage);
    path='/Users/ariannamiano/Desktop/PhD/ModelRPSPlate/strips';
    saveas(gcf,fullfile(path,['All_' num2str(t) '.svg']));
    
    elseif VIDEO == 1 %mod(t,5) ==1 && VIDEO == 1 
    rgbImage = cat(3, gridR,gridS,gridP);
    imagesc(rgbImage);
    drawnow
    frame = getframe(gcf); 
    writeVideo(writerObj, frame);
    end
end
close(writerObj)

figure
plot(gridRTotFrac,'r')
hold on
plot(gridPTotFrac,'b')
hold on
plot(gridSTotFrac,'g')


