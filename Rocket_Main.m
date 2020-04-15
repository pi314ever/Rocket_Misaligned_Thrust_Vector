%% Spin Stabilized Rocket
% Daniel Huang
% AE 140 Individual Project
% To use pre-solved presets, select the file "Rocket_Main_Preset_1.mat" when 
% prompted. 
%% Setup & Solve
clc;clear
global m l Ix Iz T G me
% Essential variables
Re = 6.378E6; % Radius of Earth [m]
G = 6.67408E-11; % Gravitational constant
rot_e = 2*pi/(24*60*60); % Earths rotation rate [rad]
init_psi = deg2rad(-120.5724); % Initial longitude
init_the = deg2rad(90-34.7420); % Initial latitude

% Check for presets
files  = dir('*.mat');
FileNames = cell(length(files),1);
for ii=1:length(files)
    FileNames{ii} = files(ii).name;
end


if ~isempty(FileNames)
    fprintf('Potential presets detected in current directory.\n')
    fprintf('Please choose from the following presets or type "0" for default values:\n')
    for ii = 1:length(FileNames)
        fprintf('\t%d) %s\n',ii,FileNames{ii})
    end
    choice = input('Please choose a file number: ');
    if choice <= length(files) && choice > 0 && choice == int64(choice)
        load(FileNames{choice})
    elseif ~choice
        fprintf('Using default values.\n')
    else
        fprintf('Invalid choice number, using in-code setup.\n')
    end
end
clear files FileNames choice % Clear unwanted variables
if (exist('G','var') && exist('m','var') && exist('me','var') && ...
        exist('l','var') && exist('Ix','var') && exist('Iz','var') &&...
        exist('T','var') && exist('al','var') && exist('tspan','var') && ...
        exist('Initials','var') && exist('Initials2','var'))
    fprintf('File has necessary variables for solution.\n')
    check = 1; % Check if solution is saved
else
    check = 0;
    fprintf('Loading in-code variables...\n')
    m = 5.49E5; % Mass of rocket [kg]
    me = 5.972E24; % Mass of earth [kg]
    r = 3.7/2; % Radius of rocket [m]
    l = 70; % Height of rocket [m]
    Ix = 0.5*m*r^2; % Ixx/Iyy of rocket about CM [kg*m^2]
    Iz = m/12*(3*r^2+l^2); % Izz of rocket about CM [kg*m^2]
    T = 7.607E6; % Thrust [N]
    al = 1E-5; % Fixed offset of thrust from nominal axis [rad]
    tspan = [0 500];
    % Initials & ODE setup
    %   INDEX   x    y    z    psi  theta phi
    %   f(t)    1    3    5    7    9     11
    %   f'(t)   2    4    6    8    10    12
    x0 = Re*sin(init_the)*sin(init_psi);
    xd0 = Re*rot_e*sin(init_the)*cos(init_psi);
    y0 = -Re*sin(init_the)*cos(init_psi);
    yd0 = Re*rot_e*sin(init_the)*sin(init_psi);
    z0 = Re*cos(init_the);
    zd0 = 0;
    psi0 = init_psi; % Precession angle [rad]
    psid0 = 0; % Precession rate [rad/s]
    the0 = init_the; % Pitch angle [rad]
    thed0 = 0; % Pitch rate [rad/s]
    phi0 = 0; % Spin angle [rad/s]
    phid0 = 0; % Spin rate [rad/s]
    Initials = [x0 xd0 y0 yd0 z0 zd0 psi0 psid0 the0 thed0 phi0 phid0]; % No spin stabilization
    Initials2 = Initials;
    Initials2(12) = pi; % Set higher spin rate [rad/s]
end

%% Check for solution
if exist('spin','var') && exist('nospin','var') && exist('al0','var')
    fprintf('Solution is pre-solved.\n')
else
    check = 0;
    % Solve 
    fprintf('Solving no-spin condition...\n')
    nospin = ode45(@(t,x)eqn(t,x,al),tspan,Initials,[]);
    fprintf('Solving spin condition...\n')
    spin = ode45(@(t,x)eqn(t,x,al),tspan,Initials2,[]);
    fprintf('Solving thrust offset angle alpha=0 without spin...\n')
    al0 = ode45(@(t,x)eqn(t,x,0),tspan,Initials,[]);
end

%% Plot 
close all
% Orientation of rocket w/o spin
orx_no = sin(nospin.y(7,:)).*sin(nospin.y(9,:));
ory_no = -sin(nospin.y(9,:)).*cos(nospin.y(7,:));
orz_no = cos(nospin.y(9,:));

% Orientation of rocket with spin
orx_sp = sin(spin.y(7,:)).*sin(spin.y(9,:));
ory_sp = -sin(spin.y(9,:)).*cos(spin.y(7,:));
orz_sp = cos(spin.y(9,:));

% Orientation of rocket with alpha = 0
orx_al0 = sin(al0.y(7,:)).*sin(al0.y(9,:));
ory_al0 = -sin(al0.y(9,:)).*cos(al0.y(7,:));
orz_al0 = cos(al0.y(9,:));

% Orientation Plots
figure
plot3(orx_no,ory_no,orz_no) % Without Spin stabilization
hold on
[a,b,c] = sphere;
mesh(a,b,c)
title('3D Orientation of Rocket Without Spin Stabilization')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
figure
plot3(orx_sp,ory_sp,orz_sp,'o') % With Spin stabilization
hold on
mesh(a,b,c)
title('3D Orientation of Rocket With Spin Stabilization')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
figure
plot3(orx_al0,ory_al0,orz_al0,'o') % Alpha = 0
hold on
mesh(a,b,c)
title('3D Orientation of Rocket With \alpha = 0')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

% Position 3D
figure
plot3(nospin.y(1,:),nospin.y(3,:),nospin.y(5,:),'b',spin.y(1,:),spin.y(3,:)...
    ,spin.y(5,:),'r',al0.y(1,:),al0.y(3,:),al0.y(5,:),'x',[x0*.9 x0],...
    [y0*.9 y0],[z0*.9 z0],':k')
title('3D Trajectory of Rocket With and Without Spin Stabilization')
legend('Without spin stabilization','With spin stabilization',...
    'With \alpha = 0','Center of Earth Reference','Location','Best')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

% Angles Comparison
figure
subplot(2,3,1) % Precession Without
plot(nospin.x,nospin.y(7,:),':r')
title('\psi Precession Angle w/o Spin')
ylabel('Angle [rad]')

subplot(2,3,4) % Precession With
plot(spin.x,spin.y(7,:),'--r',al0.x,al0.y(7,:),':o')
title('\psi Precession Angle /w Spin & \alpha = 0')
legend('With spin','With \alpha = 0','Location','Best')
xlabel('Time [s]')
ylabel('Angle [rad]')
ylim([-pi 0])

subplot(2,3,2) % Pitch Without
plot(nospin.x,nospin.y(9,:),':g')
title('\theta Pitch Angle w/o Spin')
ylabel('Angle [rad]')
ylim([0 pi])

subplot(2,3,5) % Pitch With
plot(spin.x,spin.y(9,:),'--g',al0.x,al0.y(9,:),':o') 
title('\theta Pitch Angle /w Spin & \alpha = 0')
legend('With spin','With \alpha = 0','Location','Best')
xlabel('Time [s]')
ylabel('Angle [rad]')
ylim([0 pi])

subplot(2,3,3) % Spin Without
plot(nospin.x,nospin.y(11,:),':b')
title('\phi Spin Angle w/o Spin')
ylabel('Angle [rad]')

subplot(2,3,6) % Spin With
plot(spin.x,spin.y(11,:),'--b',al0.x,al0.y(11,:),':o')
title('\phi Spin Angle /w Spin & \alpha = 0')
legend('With spin','With \alpha = 0','Location','Best')
xlabel('Time [s]')
ylabel('Angle [rad]')

% Full Screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% w3 Plots (should be constant per system)
w3ns = nospin.y(8,:).*cos(nospin.y(9,:))+nospin.y(12,:);
w3sp = spin.y(8,:).*cos(spin.y(9,:))+spin.y(12,:);
w3al0 = al0.y(8,:).*cos(al0.y(9,:))+al0.y(12,:);
figure
plot(nospin.x,w3ns,'r',spin.x,w3sp,'g',al0.x,w3al0,'ob')
title('\omega_3 Plots')
legend('No spin stabilization','Spin stabilization','Control','Location','Best')
xlabel('Time [s]')
ylabel('\omega_3 [rad]')

%% Save solutions 
if ~check
    varsav = input('Save variables? Input 0 for no, 1 for solutions only, or 2 for all variables: ');
    switch varsav
        case 1
            name = input('Please input a valid file name:\n','s');
            save(name,'nospin','spin')
            fprintf('Solutions saved in %s.mat.\n',name)
        case 2
            name = input('Please input a valid file name:\n','s');
            save(name)
            fprintf('All variables saved in %s.mat.\n',name)
        otherwise
            fprintf('No variables will be saved.\n')
    end
else
    fprintf('All variables are already saved.\n')
end
fprintf('\nProgram Complete.\n')
