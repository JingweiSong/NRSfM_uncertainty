function skelPlayData(S1,S2,Var_W)
%   Var_W [Optional]:   Drawing ellipsoid

fps = 30;



%   Modified by Jingwei
% MOCAP
pts1 = S1';
pts2 = S2';
[num_pts,num_frame] = size(pts1);
num_frame = num_frame / 3;


% Get the limits of the motion.

minY1 = min(pts1(:,1));minY2 = min(pts1(:,2));minY3 = min(pts1(:,3));
maxY1 = max(pts1(:,1));maxY2 = max(pts1(:,2));maxY3 = max(pts1(:,3));
for i = 2:num_frame
  minY1 = min([pts1(:, 3*i-2); minY1]);
  minY2 = min([pts1(:, 3*i-1); minY2]);
  minY3 = min([pts1(:, 3*i  ); minY3]);
  maxY1 = max([pts1(:, 3*i-2); maxY1]);
  maxY2 = max([pts1(:, 3*i-1); maxY2]);
  maxY3 = max([pts1(:, 3*i  ); maxY3]);
end

figure
handle(1) = plot3(pts1(:,3*1-2),pts1(:,3*1-0),pts1(:,3*1-1), '.','MarkerSize',10,'Color','b');hold on;
handle(2) = plot3(pts2(:,3*1-2),pts2(:,3*1-0),pts2(:,3*1-1), '.','MarkerSize',10,'Color','r');hold on;

if(exist('Var_W','var'))
    Var_W = 1.96*sqrt(Var_W);
    Var_W = Var_W';
    i = 1;
    for j = 1 : size(pts1,1)
        [X,Y,Z]=ellipsoid(pts1(j,3*i-2),pts1(j,3*i-0),pts1(j,3*i-1),Var_W(j,3*i-2),Var_W(j,3*i-0),Var_W(j,3*i-1));hold on;
        handle(2+j) = surf(X,Y,Z);
        handle(2+j).EdgeColor = 'none';
        handle(2+j).FaceColor = 'k';
        handle(2+j).FaceAlpha = 0.3;
        handle(2+j).FaceLighting = 'none';
    end
    legend('Prediction','Monte Carlo average')
    
    videoName = 'outvideo.avi';
    if(exist('videoName','file'))
        delete videoName.avi;
    end
    aviobj = VideoWriter(videoName); 
    aviobj.FrameRate = fps;
    open(aviobj);%Open file for writing video data  
end
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
axis on
grid on
axis([minY1,maxY1,...
        minY3,maxY3,...
        minY2,maxY2]);
view(0,45)
set(gcf, 'Position', get(0, 'Screensize'));

for i=1:num_frame
    set(handle(1), 'Xdata', pts1(:,3*i-2), 'Ydata', pts1(:,3*i-0), 'Zdata', ...
                 pts1(:,3*i-1));
    set(handle(2), 'Xdata', pts2(:,3*i-2), 'Ydata', pts2(:,3*i-0), 'Zdata', ...
                 pts2(:,3*i-1),'Color','r');
    if(exist('Var_W','var'))
        for j = 1 : size(pts1,1)
            [X,Y,Z]=ellipsoid(pts1(j,3*i-2),pts1(j,3*i-0),pts1(j,3*i-1),Var_W(j,3*i-2),Var_W(j,3*i-0),Var_W(j,3*i-1));
            set(handle(j+2), 'Xdata', X, 'Ydata', Y, 'Zdata', ...
                 Z);
        end
        frames = getframe(gcf);
        writeVideo(aviobj,frames.cdata);
    end    
    pause(1/fps);
end

if(exist('Var_W','var'))
    close(aviobj);
    disp("Video saved");
end

