% THIS IS A BODGE!!! DO NOT EXPECT TO SEE GOOD CODING PRACTICES HERE

clear
close all

%path = '~/Dropbox/TFG_MEN/Results/New/Bad linear 1 nonlinear/';
path = ''

quaternion_aux = importdata([path,'real_quaternion'],' ');
control = importdata([path,'results_control'],' ');
time = importdata([path,'results_time'],' ');

quaternion_real = zeros(size(quaternion_aux));
ver = zeros(8,3);

body_color = [.25,.25,.95];
solar_color = [.35,.35,.35];

quaternion_real(:,1) = quaternion_aux(:,4);
quaternion_real(:,2) = quaternion_aux(:,1);
quaternion_real(:,3) = quaternion_aux(:,2);
quaternion_real(:,4) = quaternion_aux(:,3);

time_length = 5;
writerObj = VideoWriter([path,'animation.avi']);
writerObj.FrameRate = length(quaternion_real)/time_length;

figure(1)
clf
open(writerObj)
for i = 1:size(quaternion_real,1)-1
    q = quaternion_real(i,:);
    x = quatrotate(quaternion_real(i,:),[1,0,0]);
    y = quatrotate(quaternion_real(i,:),[0,1,0]);
    z = quatrotate(quaternion_real(i,:),[0,0,1]);
    ver(1,:) = quatrotate(q,[1,1,1]);
    ver(2,:) = quatrotate(q,[1,1,-1]);
    ver(3,:) = quatrotate(q,[-1,1,-1]);
    ver(4,:) = quatrotate(q,[-1,1,1]);
    ver(5,:) = quatrotate(q,[1,-1,1]);
    ver(6,:) = quatrotate(q,[1,-1,-1]);
    ver(7,:) = quatrotate(q,[-1,-1,-1]);
    ver(8,:) = quatrotate(q,[-1,-1,1]);
    ver = ver*.5;
    solar1(1,:) = quatrotate(q,[.1,3,.75]);
    solar1(2,:) = quatrotate(q,[.1,3,-.75]);
    solar1(3,:) = quatrotate(q,[-.1,3,-.75]);
    solar1(4,:) = quatrotate(q,[-.1,3,.75]);
    solar1(5,:) = quatrotate(q,[.1,1,.75]);
    solar1(6,:) = quatrotate(q,[.1,1,-.75]);
    solar1(7,:) = quatrotate(q,[-.1,1,-.75]);
    solar1(8,:) = quatrotate(q,[-.1,1,.75]);
    %
    solar2(1,:) = quatrotate(q,[.1,-3,.75]);
    solar2(2,:) = quatrotate(q,[.1,-3,-.75]);
    solar2(3,:) = quatrotate(q,[-.1,-3,-.75]);
    solar2(4,:) = quatrotate(q,[-.1,-3,.75]);
    solar2(5,:) = quatrotate(q,[.1,-1,.75]);
    solar2(6,:) = quatrotate(q,[.1,-1,-.75]);
    solar2(7,:) = quatrotate(q,[-.1,-1,-.75]);
    solar2(8,:) = quatrotate(q,[-.1,-1,.75]);

    clf
    hold on
    
    % Satellite's body
    fill3([ver(1,1),ver(2,1),ver(3,1),ver(4,1)],[ver(1,2),ver(2,2),ver(3,2),ver(4,2)],[ver(1,3),ver(2,3),ver(3,3),ver(4,3)],body_color)
    fill3([ver(1,1),ver(2,1),ver(6,1),ver(5,1)],[ver(1,2),ver(2,2),ver(6,2),ver(5,2)],[ver(1,3),ver(2,3),ver(6,3),ver(5,3)],body_color)
    fill3([ver(2,1),ver(3,1),ver(7,1),ver(6,1)],[ver(2,2),ver(3,2),ver(7,2),ver(6,2)],[ver(2,3),ver(3,3),ver(7,3),ver(6,3)],body_color)
    fill3([ver(4,1),ver(3,1),ver(7,1),ver(8,1)],[ver(4,2),ver(3,2),ver(7,2),ver(8,2)],[ver(4,3),ver(3,3),ver(7,3),ver(8,3)],body_color)
    fill3([ver(5,1),ver(6,1),ver(7,1),ver(8,1)],[ver(5,2),ver(6,2),ver(7,2),ver(8,2)],[ver(5,3),ver(6,3),ver(7,3),ver(8,3)],body_color)
    fill3([ver(1,1),ver(4,1),ver(8,1),ver(5,1)],[ver(1,2),ver(4,2),ver(8,2),ver(5,2)],[ver(1,3),ver(4,3),ver(8,3),ver(5,3)],body_color)
    % Solar parel right
    fill3([solar1(1,1),solar1(2,1),solar1(3,1),solar1(4,1)],[solar1(1,2),solar1(2,2),solar1(3,2),solar1(4,2)],[solar1(1,3),solar1(2,3),solar1(3,3),solar1(4,3)],solar_color)
    fill3([solar1(1,1),solar1(2,1),solar1(6,1),solar1(5,1)],[solar1(1,2),solar1(2,2),solar1(6,2),solar1(5,2)],[solar1(1,3),solar1(2,3),solar1(6,3),solar1(5,3)],solar_color)
    fill3([solar1(2,1),solar1(3,1),solar1(7,1),solar1(6,1)],[solar1(2,2),solar1(3,2),solar1(7,2),solar1(6,2)],[solar1(2,3),solar1(3,3),solar1(7,3),solar1(6,3)],solar_color)
    fill3([solar1(4,1),solar1(3,1),solar1(7,1),solar1(8,1)],[solar1(4,2),solar1(3,2),solar1(7,2),solar1(8,2)],[solar1(4,3),solar1(3,3),solar1(7,3),solar1(8,3)],solar_color)
    fill3([solar1(5,1),solar1(6,1),solar1(7,1),solar1(8,1)],[solar1(5,2),solar1(6,2),solar1(7,2),solar1(8,2)],[solar1(5,3),solar1(6,3),solar1(7,3),solar1(8,3)],solar_color)
    fill3([solar1(1,1),solar1(4,1),solar1(8,1),solar1(5,1)],[solar1(1,2),solar1(4,2),solar1(8,2),solar1(5,2)],[solar1(1,3),solar1(4,3),solar1(8,3),solar1(5,3)],solar_color)
    % Solar panel left
    fill3([solar2(1,1),solar2(2,1),solar2(3,1),solar2(4,1)],[solar2(1,2),solar2(2,2),solar2(3,2),solar2(4,2)],[solar2(1,3),solar2(2,3),solar2(3,3),solar2(4,3)],solar_color)
    fill3([solar2(1,1),solar2(2,1),solar2(6,1),solar2(5,1)],[solar2(1,2),solar2(2,2),solar2(6,2),solar2(5,2)],[solar2(1,3),solar2(2,3),solar2(6,3),solar2(5,3)],solar_color)
    fill3([solar2(2,1),solar2(3,1),solar2(7,1),solar2(6,1)],[solar2(2,2),solar2(3,2),solar2(7,2),solar2(6,2)],[solar2(2,3),solar2(3,3),solar2(7,3),solar2(6,3)],solar_color)
    fill3([solar2(4,1),solar2(3,1),solar2(7,1),solar2(8,1)],[solar2(4,2),solar2(3,2),solar2(7,2),solar2(8,2)],[solar2(4,3),solar2(3,3),solar2(7,3),solar2(8,3)],solar_color)
    fill3([solar2(5,1),solar2(6,1),solar2(7,1),solar2(8,1)],[solar2(5,2),solar2(6,2),solar2(7,2),solar2(8,2)],[solar2(5,3),solar2(6,3),solar2(7,3),solar2(8,3)],solar_color)
    fill3([solar2(1,1),solar2(4,1),solar2(8,1),solar2(5,1)],[solar2(1,2),solar2(4,2),solar2(8,2),solar2(5,2)],[solar2(1,3),solar2(4,3),solar2(8,3),solar2(5,3)],solar_color)
    % Appendages 
    plot3([0,y(1)],[0,y(2)],[0,y(3)],'k','linewidth',5)
    plot3([0,-y(1)],[0,-y(2)],[0,-y(3)],'k','linewidth',5)
    % Body frame
    plot3([0,x(1)],[0,x(2)],[0,x(3)],'r','linewidth',1.5)
    plot3([0,y(1)],[0,y(2)],[0,y(3)],'r','linewidth',1.5)
    plot3([0,z(1)],[0,z(2)],[0,z(3)],'r','linewidth',1.5)
    % Inertial reference frame
    plot3([0,1.5],[0,0],[0,0],'g','linewidth',1.5) 
    plot3([0,0],[0,-1.5],[0,0],'g','linewidth',1.5)
    plot3([0,0],[0,0],[0,1.5],'g','linewidth',1.5)
    hold off
    xlabel('X'); ylabel('Y'), zlabel('Z')
    axis([-1,1,-1,1,-1,1]*3.5)
    axis square
    %set(gcf,'Units','centimeters')
    %set(gcf,'Position',[0,0,15,10])
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[])
    xlabel(''), ylabel(''), zlabel('')
    view([35,22.5])
    grid on
    box on
    saveas(gcf,'image','png');
    frame = im2frame(imread('image.png'));
    writeVideo(writerObj, frame);
end

for i=1:100 %Write a few more frames with the satellite at rest
    saveas(gcf,'image','png');
    frame = im2frame(imread('image.png'));
    writeVideo(writerObj, frame);
end
        
close(writerObj);
close gcf
