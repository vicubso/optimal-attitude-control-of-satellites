clear
close all
%path = '/home/vicubso/Dropbox/TFG_MEN/Results/aaa Prety bad linear!/';
%path = '/home/vicubso/Dropbox/TFG_MEN/Algorithm/satellites_nonlinear/';
path = '/home/vicubso/Dropbox/TFG_MEN/Algorithm/binaries/';
model = 1;  % 0:linear model, 1:non linear model

linear_names = {'angveloc';'control';'eulerang';'quaternion';...
    'real_angveloc';'real_eulerang';'real_quaternion'}; %;coltrol_LQR]
nonlinear_names = {'angveloc';'control';'eulerang';'quaternion';...
    'real_angveloc';'real_eulerang';'real_quaternion';'convergence';'quatnorm'};

time = importdata([path,'results_time'],' ');
angveloc = importdata([path,'results_angveloc'],' ');
control = importdata([path,'results_control'],' ');
eulerang = importdata([path,'results_eulerang'],' ');
quaternion = importdata([path,'results_quaternion'],' ');
real_angveloc = importdata([path,'real_angveloc'],' ');
real_eulerang = importdata([path,'real_eulerang'],' ');
real_quaternion = importdata([path,'real_quaternion'],' ');

error_angveloc = zeros(size(time));
error_quaternion = zeros(size(time));

for i = 1:length(time)
    error_angveloc(i) = norm(real_angveloc(i,:)-angveloc(end,:));
    error_quaternion(i) = norm(real_quaternion(i,:)-quaternion(end,:));
end

if model == 0
    control_LQR = importdata([path,'results_control_LQR'],' ');
elseif model == 1
    converg = importdata([path,'results_converg'],' ');
    quatnorm = importdata([path,'results_quatnorm'],' ');
end

pi_ticks = {'-\pi','-7\pi/8','-3\pi/4','-5\pi/8','-\pi/2','-3\pi/8',...
            '-\pi/4','-\pi/8',0,'\pi/8','\pi/4','3\pi/8','\pi/2',...
            '5\pi/8','3\pi/4','7\pi/8','\pi'};

k = 0;
% ======== ANGULAR VELOCITIES ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,angveloc(:,1),'-','linewidth',2)
xlim([time(1),time(end)])
ylim([min(min(angveloc))-pi/8,max(max(angveloc))+pi/8])
set(gca,'YTick',-pi:pi/8:pi) 
set(gca,'YTickLabel',pi_ticks)
grid on
grid minor
xlabel('t (s)')
ylabel('\omega_1 (rad/s)')

% % ======== ERROR ANG. VEL. ========
% k = k+1;
% figure(k)
% clf
% set(gcf,'Units','centimeters')
% set(gcf,'Position',[0,0,15,10])
% set(gca, 'FontName', 'Arial')
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf, 'resize', 'off')
% hold on
% plot(time,error_angveloc,'-','linewidth',2)
% plot([0,time(end)],[0,0],'k--','linewidth',1)
% hold off
% xlim([time(1),time(end)])
% ylim([-pi/8,max(error_angveloc)+pi/8])
% set(gca,'YTick',-pi:pi/8:pi) 
% set(gca,'YTickLabel',pi_ticks)
% grid on
% grid minor
% xlabel('t (s)')
% ylabel('|\omega-\omega_T|(rad/s)')

% ========= CONTROL ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
hold on
plot(time(1:end-1),control(:,1),'-','linewidth',2)
plot([0,time(end)],[1,1],'k--','linewidth',1.5)
plot([0,time(end)],[-1,-1],'k--','linewidth',1.5)
hold off
xlim([time(1),time(end-1)])
ylim([-1.25,1.25])
grid on
grid minor
xlabel('t (s)')
ylabel('u_1 (Nm)')

% ======== EULER ANGLES ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,eulerang(:,3),'-','linewidth',2)
xlim([time(1),time(end)])
xlim([time(1),time(end)])
ylim([min(min(eulerang))-pi/8,max(max(eulerang))+pi/8])
set(gca,'YTick',-pi:pi/8:pi) 
set(gca,'YTickLabel',pi_ticks)
grid on
grid minor
xlabel('t (s)')
ylabel('\psi (rad)')


% ======== QUATERNION ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,quaternion,'-','linewidth',2)
xlim([time(1),time(end)])
ylim([min(min(quaternion))-0.2,max(max(quaternion))+0.5])
grid on
grid minor
xlabel('t (s)')
ylabel('q')
legend('\epsilon_1','\epsilon_1','\epsilon_3','\eta')

% % ======== ERROR QUATERNION ========
% k = k+1;
% figure(k)
% clf
% set(gcf,'Units','centimeters')
% set(gcf,'Position',[0,0,15,10])
% set(gca, 'FontName', 'Arial')
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf, 'resize', 'off')
% hold on
% plot(time,error_quaternion,'-','linewidth',2)
% plot([0,time(end)],[0,0],'k--','linewidth',1)
% hold off
% xlim([time(1),time(end)])
% ylim([-0.2,max(error_quaternion)+0.2])
% grid on
% grid minor
% xlabel('t (s)')
% ylabel('|q-q_T|')

% ======== REAL ANGULAR VELOCITIES ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,real_angveloc(:,1),'-','linewidth',2)
xlim([time(1),time(end)])
ylim([min(min(real_angveloc))-pi/8,max(max(real_angveloc))+pi/8])
set(gca,'YTick',-pi:pi/8:pi) 
set(gca,'YTickLabel',pi_ticks)
grid on
grid minor
xlabel('t (s)')
ylabel('\omega_1 (rad/s)')

% ======== REAL EULER ANGLES ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,real_eulerang(:,3),'-','linewidth',2)
xlim([time(1),time(end)])
xlim([time(1),time(end)])
ylim([min(min(real_eulerang))-pi/8,max(max(real_eulerang))+pi/8])
set(gca,'YTick',-pi:pi/8:pi) 
set(gca,'YTickLabel',pi_ticks)
grid on
grid minor
xlabel('t (s)')
ylabel('\psi (rad)')

% ======== REAL QUATERNION ========
k = k+1;
figure(k)
clf
set(gcf,'Units','centimeters')
set(gcf,'Position',[0,0,15,10])
set(gca, 'FontName', 'Arial')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'resize', 'off')
plot(time,real_quaternion,'-','linewidth',2)
xlim([time(1),time(end)])
ylim([min(min(real_quaternion))-0.2,max(max(real_quaternion))+0.5])
grid on
grid minor
xlabel('t (s)')
ylabel('q')
legend('\epsilon_1','\epsilon_1','\epsilon_3','\eta')

if model == 0
%     % ======== CONTROL FROM LQR ========
%     k = k+1;
%     figure(k)
%     clf
%     plot(time,control_LQR,'-','linewidth',2)
%     xlim([time(1),time(end-1)])
%     grid on
%     xlabel('t (s)')
%     ylabel('u (Nm)')
%     legend('u_1','u_2','u_3')
elseif model == 1
    % ======== CONVERGENCE =======
    k = k+1; %TODO log axis x
    figure(k)
    clf
    set(gcf,'Units','centimeters')
    set(gcf,'Position',[0,0,15,10])
    set(gca, 'FontName', 'Arial')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf, 'resize', 'off')
    loglog(converg,'k-','linewidth',2)
    xlim([1,length(converg)])
    ylim([1e-4,max(max(converg))+0.2])
    grid on
    grid minor
    xlabel('Iterations')
    ylabel('max\{v_k R_k\}')
    
    % ======== NORM OF QUATERNION
    k = k+1;
    figure(k)
    clf
    set(gcf,'Units','centimeters')
    set(gcf,'Position',[0,0,15,10])
    set(gca, 'FontName', 'Arial')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf, 'resize', 'off')
    plot(time,quatnorm,'k-','linewidth',2)
    xlim([time(1),time(end)])
    ylim([min(min(quatnorm))-2e-6,max(max(quatnorm))+2e-6])
    grid on
    grid minor
    xlabel('t (s)')
    ylabel('|q|')
end

flag = input(' ');
if flag == 1
    if model == 0
        names = linear_names;
    elseif model == 1
        names = nonlinear_names;
    end
    for i = 1:length(names)
        saveas(figure(i),[path,'figures_',names{i}],'fig')
        saveas(figure(i),[path,'figures_',names{i}],'png')
    end
    close all
else
    close all
end