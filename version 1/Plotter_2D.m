%% makes a heatmap of one frame

dimX = 100;
dimY = 100;

X = zeros(dimX,dimY);

for i = 1:dimX
    for j = 1:dimY
        X(i,j) = grandi_2D_3c7_s1s2_1Hz_def_100deltaT(1100,(i-1)*dimX + j);
    end    
end

HeatMap(X, 'Symmetric', 'false');

%% creates a matlab movie

numFrames = 1000;
dimX = 100;
dimY = 100;

X = zeros(dimX,dimY);
A = moviein(numFrames);
set(gca, 'NextPlot','replacechildren');
axis([0 dimX 0 dimY -80 10 -80 10]);
axis off;

for h = 1:(numFrames)
    for i = 1:dimX
        for j = 1:dimY
            X(i,j) = grandi_2D_3c7_s1s2_1Hz_def_453deltaT_BigS2_400x400(h,(i-1)*dimX + j);
        end    
    end
    
    surf(X,'LineStyle','none');
    shading(gca,'interp');
    view(0,90); %this makes the view from the top, if you comment it out you get the 3d plot
    A(h) = getframe;

end

%movie(A,1,5);

%% creates AVI movie

writerObj = VideoWriter('grandi_2D_3c7_s1s2_1Hz_def_453deltaT_bigS2_400x400.avi');
writerObj.FrameRate = 20;
open(writerObj);

numFrames = 1000;
dimX = 100;
dimY = 100;

X = zeros(dimX,dimY);
set(gca, 'NextPlot','replacechildren');
set(gcf,'Renderer','zbuffer');
axis([0 dimX 0 dimY -80 30 -80 30]);
axis off;

for h = 1:(numFrames)
    for i = 1:dimX
        for j = 1:dimY
            X(i,j) = grandi_2D_3c7_s1s2_1Hz_def_453deltaT_BigS2_400x400(h,(i-1)*dimX + j);
        end    
    end
    
    surf(X,'LineStyle','none');
    view(0,90);
    A = getframe;
    writeVideo(writerObj, A);
    
end

close(writerObj);

