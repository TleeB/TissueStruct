%# figure
figure, set(gcf, 'Color','white')

Z = peaks; surf(Z);  axis tight
set(gca, 'nextplot','replacechildren', 'Visible','off');

%# preallocate
nFrames = 200;
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);

%# create movie
for k=1:nFrames
   surf(sin(2*pi*k/20)*Z, Z)
   mov(k) = getframe(gca);
end
close(gcf)

%# save as AVI file, and open it using system video player
movie2avi(mov, 'myPeaks1.avi', 'compression','None', 'fps',15);
%winopen('myPeaks1.avi')

N = reshape(testGrandi(1,:), 5,4);

