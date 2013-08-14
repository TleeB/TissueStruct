% Simple test of the Matlab's ability to make a movie without putting
% anything on the screen.  (Prevents work overalapping with the
% movie-plotting windwow from appearing in the movie.)
filename = 'movie_making_test.avi';
mov = avifile(filename);
fig = figure('visible','off');
% The following determines the figure size.  To obtain these numbers,
% set 'off' to 'on' temporarily in the line above, the run the code
% with a small number of frames.  %Expand the figure as desired, then
% execute the matlab command, 'get(fig)'. The four values printed out for
% the 'Position' command are the ones to use here.  Then set the 'on' back
% to 'off' in the above command.
set(fig,'Position',[60,0,643,591]);
clf;
for fr = 1:60
    x = ones(128,1)*(-64:63);
    y = (-64:63)'*ones(1,128);
    z = exp(-(x.*x+y.*y)/(2*fr*fr));
    imagesc(z); colorbar;
    xlabel('x'); ylabel('y'); title(sprintf('Gaussian for frame no. %i',fr));
    drawnow;
    mov = addframe(mov,fig);
    fprintf('Plotted frame no. %i\n',fr);
end
mov = close(mov);