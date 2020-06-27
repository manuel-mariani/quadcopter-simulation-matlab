%% Mathematical Modeling and Simulation of Quadcopter Drone
%  Manuel Mariani, Daniele Manes; Ancona 2018
clc
clear
close all

%% Variables & Parameters Declaration
% Parameters
Deltat = 0.01;       % Value of discretization of the time interval (s)
m = 0.45;            % Mass of drone (Kg)
l = 0.23;            % Length of drone arms, from the center (m)
b = 7.5e-7 ;         % Drag constant
L = 3e-6;            % Lift constant
A = 0.25/m*eye(3);   % Aerodynamical effects matrix
cost = l*L;
% Inertia matrix
Ixx = 5e-3;
Iyy = 5e-3;
Izz = 8e-3;
I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
Ir = 6e-5;          


%    Variables
% w is the matrix of angular velocities of the 4 motors x Time interval
% (rad/s)
w = zeros(4,1);
%   Initial Conditions, supposedly zero; Preallocation
xi = zeros(3, 2);   % Relative position (x,y,z) of drone, relative to Earth Frame
xii = zeros(3,2);   % Derivate of xi
eta = zeros(3, 2);  % Attitudes of drone (roll, pitch, yaw)
tau = zeros(3,1);   % Torques of drone, relative to body frame attitudes
T = 0;              % Thrust of drone on the z axis
vel = zeros(3, 2);  % Angular velocities relative to the Body Frame
Winv = zeros(3,3);  % Tranformation matrix (inverted)
R = zeros(3,1);     % Rotation vector of drone's z axis
Rx = zeros(3,3,2);  % x Rotation matrix
Ry = zeros(3,3,2);  % y Rotation matrix
Rz = zeros(3,3,2);  % z Rotation matrix

%% w custom values
% 613 for almost stationary flight
wBase = 613;        % Base value for motor's speeds
deltaW = 0.5;       % Step size for controlling the motors

w(1,1) =  wBase;
w(2,1) =  wBase;
w(3,1) =  wBase;
w(4,1) =  wBase;

%% Plot preparation
k = 0;              % Counter
Simulate = true;    % when false, interrupts simulation
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
set(f1,'WindowKeyPressFcn',@KeyDown, 'WindowKeyReleaseFcn', @KeyUp); 
% Drone shape definition
Drone = [ 1 -1 0 0 0        
          0 0 0 1 -1        
          0 0 0 0 0];       
RotatedDrone = Drone;
% Main plot for 3D simulation
droneView = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
p = plot3(RotatedDrone(1,:,1),RotatedDrone(2,:,1),RotatedDrone(3,:,1),'k.-');
animLine = animatedline('MaximumNumPoints',1000,...
                        'Color','r');
xlabel('x');ylabel('y');zlabel('z')
axisLim = 30;
axis([-axisLim axisLim -axisLim axisLim -axisLim axisLim])
grid on;
grid minor;

% Annotations for displaying numerical values
annotMotors = annotation('textbox',[3/4 3/4 1/4-0.001 1/4], ...
                         'FontSize', 15, ...
                         'BackgroundColor', [1 1 1], ...
                         'LineWidth', 0.7, ...
                         'Margin', 10,...
                         'FontUnits','normalized');
annotAngles = annotation('textbox',[3/4 1/4 1/4-0.001 2/4], ...
                         'FontSize', 15, ...
                         'BackgroundColor', [1 1 1], ...
                         'LineWidth', 0.7, ...
                         'Margin', 10, ...
                         'FontUnits','normalized');
annotControls = annotation('textbox',[3/4 0 1/4 1/4], ...
                         'FontSize', 15, ...
                         'BackgroundColor', [1 1 1], ...
                         'LineWidth', 0.7, ...
                         'Margin', 10,...
                         'FontUnits','normalized');
controlsGuideTxt = ['Controls:' newline ...
                        '  W and S :  Pitch' newline ...
                        '  A and D :  Roll' newline ...
                        '  Left and Right: Yaw' newline ...
                        '  Up and Down : All motors'' speed' newline ...
                        '  R : Reset simulation' newline ...
                        '  C : Close simulation'];
set(annotControls,'String', controlsGuideTxt);

% UIcontrol
annotSliders = annotation('textbox',[0 0 3/4 1/4], ...
                         'FontSize',15, ...
                         'BackgroundColor', [1 1 1], ...
                         'FontUnits','normalized',...
                         'String', ...
                         [newline 'Zoom:' newline newline 'Base speed:' ...
                         newline newline 'Speed increment:']);
sliderZoom = uicontrol('Style','slider', 'Min',0.1, 'Value',0.5,...
                         'Units','normalized',...
                         'Position', [0.12 0.19 0.2 0.025]);
sliderSpeed = uicontrol('Style','slider', 'Min',0, 'Value',613,...
                         'Units','normalized',...
                         'Position', [0.12 0.125 0.2 0.025], ...
                         'Max',1000, 'Callback',@cbkSliderSpeed);
sliderIncrement = uicontrol('Style','slider', 'Min',0, 'Value',0.5,...
                         'Units','normalized',...
                         'Position', [0.12 0.06 0.2 0.025], ...
                         'Max',1, 'Callback',@cbkSliderIncrement);
buttonRecord = uicontrol('String', 'Record Video',...
                         'Units','normalized',...
                         'Position', [0.65 0.19 0.06 0.025],...
                         'Callback', @cbkButtonRecord);

recordVideo = false;


%% Simulation & Numeric implementation
while Simulate
    k = k+1;
    
    phi = eta(1,1);
    theta = eta(2,1);
    psi = eta(3,1);
    
    w1square = w(1)^2;
    w2square = w(2)^2;
    w3square = w(3)^2;
    w4square = w(4)^2;
    % Angular momentum
    tau(:) = [cost*(w4square-w2square);
              cost*(w3square-w1square);
              b*(w1square - w2square + w3square - w4square)];
    % Thrust on drone's z axis
    T = L*(w1square + w2square + w3square + w4square);
    % Angular velocities
    Gamma = Ir * cross(vel(:,1), [0;0;1]) * (w(1,1) - w(2,1) + w(3,1) - w(4,1));
    vel(:,1+1) = Deltat*( I\(-cross( vel(:,1), I*vel(:,1)) -Gamma +tau(:)) ) + vel(:,1);
    % Transformation matrix from Body to Earth frame
    Winv(:,:) = [ 1 sin(phi)*tan(theta) cos(phi)*tan(theta);
                  0 cos(phi) -sin(phi);
                  0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
    % Attitudes
    eta(:,1+1) = Deltat*( Winv(:,:)*vel(:,1)) + eta(:,1);
    %                   Rotation matrices
    % Defining 3rd column of rotation matrix
    R(:) = [cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
            sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
            cos(theta)*cos(phi)];
    % Defining the three rotation matrices for 3d plotting
    Rx(:,:,1) = [ 1 0 0;
                  0 cos(phi) sin(phi);
                  0 -sin(phi) cos(phi)];
    Ry(:,:,1) = [cos(theta) 0 -sin(theta);
                 0 1 0;
                 sin(theta) 0 cos(theta)];
    Rz(:,:,1) = [cos(psi) sin(psi) 0;
                 -sin(psi) cos(psi) 0;
                 0 0 1];
  
    % Positions (absolute)
    xii(:,1+1) = Deltat*( -[0;0;10] + T/m * R(:) -A*xii(:,1)) + xii(:,1);
    xi(:,1+1) = (Deltat*xii(:,1)) + xi(:,1);  
    
    %% Plotting
    % 3d drone plot
    if mod(k,2) == 0
       axisLim = 10/sliderZoom.Value;
       pause(0.001)
       axis([-axisLim+xi(1,1) axisLim+xi(1,1) -axisLim+xi(2,1)...
              axisLim+xi(2,1) -axisLim+xi(3,1) axisLim+xi(3,1)]);
        
        
        % Rotation of drone
        RotatedDrone(:,:) = Rx(:,:,1)'*Ry(:,:,1)'*Rz(:,:,1)*(Drone) + xi(:,1);
        set(p, 'XData', RotatedDrone(1,:), ...
               'YData', RotatedDrone(2,:), ...
               'ZData', RotatedDrone(3,:));          
        % Drone's trail
        addpoints(animLine,xi(1,1),xi(2,1),xi(3,1));
        drawnow;
        pause(0.001);
      
    end
    
    % Update text annotations
    if mod(k,20) == 0
      txtMotors = ['Motors'' speeds (rad/s)' newline ...
          '    \omega_1 =' num2str(w(1)) newline ... 
          '    \omega_2 =' num2str(w(2)) newline ...
          '    \omega_3 =' num2str(w(3)) newline ...
          '    \omega_4 =' num2str(w(4)) newline];
      txtAngles = ['Angles'' values (rad)' newline ...
          '    \phi =' num2str(round(phi,3)) newline ...
          '    \theta =' num2str(round(theta,3)) newline ...
          '    \psi =' num2str(round(psi,3)) newline ...
          'Torques'' values' newline ...
          '    \tau_{\phi} =' num2str(round(tau(1),3)) newline ...
          '    \tau_{\theta} =' num2str(round(tau(2),3)) newline ...
          '    \tau_{\psi} =' num2str(round(tau(3),3)) newline ...
          'Angular velocities' newline ...
          '    \nu_{\phi} =' num2str(round(vel(1,1),3)) newline ...
          '    \nu_{\theta} =' num2str(round(vel(2,1),3)) newline ...
          '    \nu_{\psi} =' num2str(round(vel(3,1),3)) newline];
      % Set new plot values
      set(annotMotors,'String',txtMotors);
      set(annotAngles,'String',txtAngles);
      drawnow limitrate
    end
    if mod(k,5) == 0
      if recordVideo
        writeVideo(VW,getframe(f1));
      end
    end
    
    %% Setting values of next iteration
    vel(:,1) = vel(:,2);
    eta(:,1) = eta(:,2);
    xii(:,1) = xii(:,2);
    xi(:,1) = xi(:,2);
end

%% Key functions
function KeyDown(~,key)
    switch key.Key
        case 'w'
            evalin('base','w(1) = w(1) -deltaW;');
            evalin('base','w(3) = w(3) +deltaW;');
        case 's'
            evalin('base','w(1) = w(1) +deltaW;');
            evalin('base','w(3) = w(3) -deltaW;');
        case 'a'
            evalin('base','w(2) = w(2) +deltaW;');
            evalin('base','w(4) = w(4) -deltaW;');   
        case 'd'
            evalin('base','w(2) = w(2) -deltaW;');
            evalin('base','w(4) = w(4) +deltaW;');
        case 'uparrow'
            evalin('base','w = w +deltaW*10;');
        case 'downarrow'
            evalin('base','w = w -deltaW*10;');
        case 'leftarrow'
            evalin('base', 'w = [w(1)-deltaW;w(2)+deltaW;w(3)-deltaW;w(4)+deltaW]');
        case 'rightarrow'
            evalin('base', 'w = [w(1)+deltaW;w(2)-deltaW;w(3)+deltaW;w(4)-deltaW]');  
        case 'r'
            evalin('base', 'wBase = 613;');
            evalin('base', 'w = [wBase;wBase;wBase;wBase];');
            evalin('base', 'xi = zeros(3, 2); ');
            evalin('base', 'xii = zeros(3,2); ');
            evalin('base', 'eta = zeros(3, 2);');
            evalin('base', 'tau = zeros(3,1); ');
            evalin('base', 'T = 0;            ');
            evalin('base', 'vel = zeros(3, 2);');
            evalin('base', 'Winv = zeros(3,3);');
            evalin('base', 'R = zeros(3,1);   ');
            evalin('base', 'Rx = zeros(3,3,2);');
            evalin('base', 'Ry = zeros(3,3,2);');
            evalin('base', 'Rz = zeros(3,3,2);');
            evalin('base', 'clearpoints(animLine)');
        case 'c'
            evalin('base', 'Simulate = false;');
            evalin('base', 'close');
    end
end
function KeyUp(~,~)
    evalin('base', 'w = [wBase;wBase;wBase;wBase];');
end
function cbkSliderSpeed(~,~)
    evalin('base', 'wBase = sliderSpeed.Value;');
    evalin('base', 'w = [wBase;wBase;wBase;wBase];');
end
function cbkSliderIncrement(~,~)
    evalin('base', 'deltaW = sliderIncrement.Value');
end
function cbkButtonRecord(~,~)
    evalin('base', ['if ~recordVideo VW = VideoWriter(''drone_capture.mp4'','...
                    '''MPEG-4''); open(VW); recordVideo = true;'...
                    'buttonRecord.BackgroundColor = [1 0.4 0.4];' ...
                    'buttonRecord.String = ''Stop Recording'';' ...
                    'else recordVideo = false; close(VW);' ...
                    'buttonRecord.String = ''Record Video'';' ...
                    'buttonRecord.BackgroundColor = [.94 .94 .94];' ...
                    'end']);
end