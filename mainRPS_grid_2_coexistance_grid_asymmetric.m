
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
Ntrial = 100]\

%set density 0 for lowest, 1 for medium and 2 for high;
density=1;




 if density ==0  %4 pixels of distance between consecutive colonies (4500 um)
        
     %populate R
         count(1) = 1;
         count2(1) =1;
        for n = 1:round(size(gridR,1)/15)  %populate the first column
            for m = 1:round(size(gridR,1)/15)
            ProbR = 0.5;
             x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1; 
              else
               gridP(count(n),count2(m)) = 1;   
              end  
              count(n+1) = count(n) +15;
              count2(m+1) = count2(m) +15;
            end
        end
        
        
         count(1) = 11;
         count2(1) =6;
        for n = 1:round(size(gridR,1)/15)  %populate the second column
            for m = 1:round(size(gridR,1)/15)
             x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1;
                else
               gridP(count(n),count2(m)) = 1; 
              end  
              count(n+1) = count(n) +15;
              count2(m+1) = count2(m) +15;
            end
        end
    
        
         count(1) = 6;
         count2(1) =11;
        for n = 1:round(size(gridS,1)/15)  %populate the third column
            for m = 1:round(size(gridS,1)/15)
             x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1;  
                else
               gridP(count(n),count2(m)) = 1; 
              end  
              count(n+1) = count(n) +15;
              count2(m+1) = count2(m) +15;
            end
        end

         %populate P
         count(1) = 6;
         count2(1) =1;
        for n = 1:round(size(gridP,1)/15)  %populate the first column
            for m = 1:round(size(gridP,1)/15)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end
        
        
         count(1) = 1;
         count2(1) =6;
        for n = 1:round(size(gridP,1)/15)  %populate the first column
            for m = 1:round(size(gridP,1)/15)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end
    
        
         count(1) = 11;
         count2(1) =11;
        for n = 1:round(size(gridP,1)/15)  %populate the first column
            for m = 1:round(size(gridP,1)/15)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end
        
         %populate S
         count(1) = 11;
         count2(1) =1;
        for n = 1:round(size(gridR,1)/15)  %populate the first column
            for m = 1:round(size(gridR,1)/15)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end
        
        
         count(1) = 6;
         count2(1) =6;
        for n = 1:round(size(gridR,1)/15)  %populate the first column
            for m = 1:round(size(gridR,1)/15)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end
    
        
         count(1) = 1;
         count2(1) =11;
        for n = 1:round(size(gridR,1)/15)  %populate the first column
            for m = 1:round(size(gridR,1)/15)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +15;
             count2(m+1) = count2(m) +15;
            end
        end

 end



if density ==1  %2 pixels of distance between consecutive colonies (2250 um)
       %populate S
         count(1) = 1;
         count2(1) =1;
        for n = 1:round(size(gridS,1)/9)  %populate the first column
            for m = 1:round(size(gridS,1)/9)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +9;
             count2(m+1) = count2(m) +9;
            end
        end
        
        
         count(1) = 7;
         count2(1) =4;
        for n = 1:round(size(gridS,1)/9)  %populate the first column
            for m = 1:round(size(gridS,1)/9)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +9;
             count2(m+1) = count2(m) +9;
            end
        end
    
        
         count(1) = 4;
         count2(1) =7;
        for n = 1:round(size(gridS,1)/9)  %populate the first column
            for m = 1:round(size(gridS,1)/9)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +9;
             count2(m+1) = count2(m) +9;
            end
        end

         %populate P
         count(1) = 4;
         count2(1) =1;
            for n = 1:round(size(gridP,1)/9)  %populate the first column
                for m = 1:round(size(gridP,1)/9)
                 gridP(count(n),count2(m)) = 1;  
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end
        
        
         count(1) = 1;
         count2(1) =4;
        for n = 1:round(size(gridP,1)/9)  %populate the first column
                for m = 1:round(size(gridP,1)/9)
                 gridP(count(n),count2(m)) = 1;  
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end
    
        
         count(1) = 7;
         count2(1) =7;
        for n = 1:round(size(gridP,1)/9)  %populate the first column
                for m = 1:round(size(gridP,1)/9)
                 gridP(count(n),count2(m)) = 1;  
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end
        
         %populate R
         count(1) = 7;
         count2(1) =1;
          for n = 1:round(size(gridR,1)/9)  %populate the first column
                for m = 1:round(size(gridR,1)/9)
                  ProbR = 0.2;
             x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1; 
              else
               gridP(count(n),count2(m)) = 1;   
              end    
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end
        
        
         count(1) = 4;
         count2(1) =4;
        for n = 1:round(size(gridR,1)/9)  %populate the first column
                for m = 1:round(size(gridR,1)/9)
                 x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1; 
              else
               gridP(count(n),count2(m)) = 1;   
              end      
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end
        
         count(1) = 1;
         count2(1) =7;
        for n = 1:round(size(gridR,1)/9)  %populate the first column
                for m = 1:round(size(gridR,1)/9)
                 x = rand; %rand variable between 0 and 1
              if x < ProbR
               gridR(count(n),count2(m)) = 1; 
              else
               gridP(count(n),count2(m)) = 1;   
              end      
                 count(n+1) = count(n) +9;
                 count2(m+1) = count2(m) +9;
                end
            end

end

if density ==2  %1 pixels of distance between consecutive colonies (1125 um)
      %populate S
         count(1) = 1;
         count2(1) =1;
        for n = 1:round(size(gridS,1)/6)  %populate the first column
            for m = 1:round(size(gridS,1)/6)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
        
         count(1) = 5;
         count2(1) =3;
        for n = 1:round(size(gridS,1)/6)  %populate the first column
            for m = 1:round(size(gridS,1)/6)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
    
        
         count(1) = 3;
         count2(1) =5;
        for n = 1:round(size(gridS,1)/6)  %populate the first column
            for m = 1:round(size(gridS,1)/6)
             gridS(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end

         %populate P
         count(1) = 3;
         count2(1) =1;
        for n = 1:round(size(gridP,1)/6)  %populate the first column
            for m = 1:round(size(gridP,1)/6)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
        
         count(1) = 1;
         count2(1) =3;
       for n = 1:round(size(gridP,1)/6)  %populate the first column
            for m = 1:round(size(gridP,1)/6)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
    
        
         count(1) = 5;
         count2(1) =5;
        for n = 1:round(size(gridP,1)/6)  %populate the first column
            for m = 1:round(size(gridP,1)/6)
             gridP(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
        
         %populate R
         count(1) = 5;
         count2(1) =1;
        for n = 1:round(size(gridR,1)/6)  %populate the first column
            for m = 1:round(size(gridR,1)/6)
             gridR(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
        
        
         count(1) = 3;
         count2(1) =3;
        for n = 1:round(size(gridR,1)/6)  %populate the first column
            for m = 1:round(size(gridR,1)/6)
             gridR(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
        
         count(1) = 1;
         count2(1) =5;
       for n = 1:round(size(gridR,1)/6)  %populate the first column
            for m = 1:round(size(gridR,1)/6)
             gridR(count(n),count2(m)) = 1;  
             count(n+1) = count(n) +6;
             count2(m+1) = count2(m) +6;
            end
        end
end



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
InitRed=sum(sum(gridR));
InitBlue=sum(sum(gridP));
InitGreen=sum(sum(gridS));
tot = InitRed +InitBlue+InitGreen;
disp(InitRed)
disp(InitBlue)
disp(InitGreen)
%R is red
%P is blue
%S is green

%define probability of killing, they are the max value when the cell is
%fully surrounded
%lower value correspond to the probability of dying if no killers are
%around
%top value is the highest probability if all cells are around

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

%pP =  linspace(0.05,0.5,8);
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
     if mod(t,5) ==1 && PLOT == 1
    rgbImage = cat(3, gridR,gridS,gridP);  %need to get it as a 3D matrix to get them to overlap
    imagesc(rgbImage);
    title(['time:',num2str(t),'Red=',num2str(round(gridRTotFrac(t),2)),'%, Blue=',num2str(round(gridPTotFrac(t),2)),'%, Green=',num2str(round(gridSTotFrac(t),2)),'%'])
     drawnow
    pause(0.2)
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


