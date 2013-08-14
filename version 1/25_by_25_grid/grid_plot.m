grid_values = importdata('grid_file.dat');           % A 5-by-5 matrix of random values from 0 to 1
label_flag = 1;
figure(1);
imagesc(sign(grid_values));            % Create a colored plot of the matrix values
map = hsv(20);
map2 = spring(25);
colormap([map(15,:);map2(12,:)]);  % Change the colormap to gray (so higher values are
                         %   black and lower values are white
 %set(gca,'YDir','normal')                        

if(label_flag)
textStrings = num2str(grid_values(:),'%d');  % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x,y] = meshgrid(1:length(grid_values));   % Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      % Plot the strings
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  % Get the middle value of the color range
textColors = repmat(grid_values(:) < midValue,1,3);  % Choose white or black for the
                                             %   text color of the strings so
                                             %  they can be easily seen over
                                             %  the background color
set(hStrings,{'Color'},num2cell(textColors,2));  % Change the text colors
end

%view(-37.5,30);
[nX,nY,~] = size (grid_values);
xlim([0.5,nX+.5])
ylim([0.5,nY+.5])

set(gca,'xtick',0.5:nX+.5)
set(gca,'ytick',0.5:nY+.5)
set(gca,'XTickLabel','', 'YTickLabel','')
grid
axis square

%%
myo = importdata('testGrandi.dat');
fib = importdata('testGrandi2.dat');

[rows, cols]=size(grid_values);
full_file=zeros(length(myo),1);
for R = 1:rows
    for C= 1:cols
        if grid_values(R,C) <0
            full_file = horzcat(full_file, fib(:,abs(grid_values(R,C))));
        else
            full_file = horzcat(full_file, myo(:,abs(grid_values(R,C))));
        end
    end
    
end

full_file=full_file(:,2:end);%remove first column of zeros
%%
figure(2)
%# figure
%figure, set(gcf, 'Color','white')
%axis manual

%set(gca, 'nextplot','replacechildren', 'Visible','off');
set(gca, 'NextPlot','replacechildren', 'YDir', 'reverse');%change the y axis direction

axis([1 cols 1 rows -90 50 -90 50]);
axis square
%hold on
%# preallocate
nFrames = length(myo);
mov(1:10:nFrames) = struct('cdata',[], 'colormap',[]);

%# create movie
for k=1:10:nFrames
    k
   %imagesc(reshape(full_file(k,:), cols,rows).')
   %mov(k) = getframe(gca);
   surf(reshape(full_file(k,:), cols,rows)');%,'LineStyle','none');
    shading(gca,'interp');
    view(0, 90)
    %view(15,30); %this makes the view from the top, if you comment it out you get the 3d plot
    mov(k) = getframe;

end
%close(gcf)
%movie(mov, -1)

%# save as AVI file, and open it using system video player
%movie2avi(mov, 'testgrandi.avi', 'compression','None', 'fps',100);
%winopen('myPeaks1.avi')
%% creates AVI movie
figure(3)
writerObj = VideoWriter('finally_propagation.avi');
writerObj.FrameRate = 20;
open(writerObj);

set(gca, 'NextPlot','replacechildren');
set(gcf,'Renderer','zbuffer');
axis([1 cols 1 rows -90 50 -90 50]);
axis square
%axis off;

nFrames = length(myo);
mov(1:10:nFrames) = struct('cdata',[], 'colormap',[]);
for k=1:10:nFrames
    surf(reshape(full_file(k,:), cols,rows)','LineStyle','none');
    shading(gca,'interp');
    %view(37.5*2,30);
    A = getframe;
    writeVideo(writerObj, A);
    
end

close(writerObj);

