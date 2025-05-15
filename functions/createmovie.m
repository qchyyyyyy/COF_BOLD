function createmovie(F,name,framerate)
v = VideoWriter(name,'MPEG-4');
v.Quality = 100;
%v.LosslessCompression = true;
v.FrameRate = framerate;% frames per second
open(v);
writeVideo(v,F);
close(v);
end