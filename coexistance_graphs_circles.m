close all
array = size(steadystate,1);
trials = size(steadystate,2);

%steadystate(7:9,:,:) = steadystate2(:,:,:);

%red,blue,green
winnerRed = zeros(array,trials);
winnerBlue = zeros(array,trials);
winnerGreen = zeros(array,trials);
coexistance = zeros(array,trials);

winR=zeros(array,1);
winB=zeros(array,1);
winG=zeros(array,1);
coex=zeros(array,1);


for a=1:array
    for t=1:trials
        
        if steadystate(a,t,1)>0.9
            winnerRed(a,t) =1;
            winR(a) =winR(a)+1;
        elseif steadystate(a,t,2)>0.9
            winnerBlue(a,t) =1;
            winB(a) =winB(a)+1;
        elseif steadystate(a,t,3)>0.9
            winnerGreen(a,t) =1;
            winG(a) =winG(a)+1;
        else
            coexistance(a,t) =1;
            coex(a) = coex(a)+1;
        end
        
    end
end


together =[winG';winB';winR';coex'];

x=0.1:0.05:0.5;
figure
bar(x,together','stacked')
xlabel('Probability of death of Strain P (Blue)')
legend('Green wins', 'Blue wins','Red wins','coexistance')
title('pR = 0.1, pS=0.5')


figure
rgbImage = cat(3, winnerRed, winnerGreen,winnerBlue);  %need to get it as a 3D matrix to get them to overlap
imagesc(rgbImage);



