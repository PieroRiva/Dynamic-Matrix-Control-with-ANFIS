%% Blimp control by ANFIS + MPC
% By Piero Alessandro Riva Riquelme, 24yo. 
% Undergraduate student. 
% Universidad de Concepción, Chile.
% pieroriva@udec.cl

% Blimp equations
% Based on work of Michale Frye, Stephen Gammon and Chunjiang Qian
% "The 6-DOF dynamic model and simulation of the tri-turbofan remote-controlled airship".
% Based on work of S.B.V. Gomes 
% "An investigation of the flight dynamics of airships with application to the YEZ-2A".
% "Airship dynamic modeling for autonomous operation".
%% Step 0: Cleaning
clear all
close all
clc
fprintf('\nANFIS SIMULATION STARTED...');
tic

%% Step 1: Process definition (6DOF blimp)
fprintf('\n\nStep 1: Process definition (6DOF blimp)...');
	% a) Physical constants
    rhoAir      = 1.205; % Density of air at NTP (20°C, 1atm)
    rhoHe       = 0.1664; % Density of helium at NTP (20°C, 1atm)
    gacc        = 9.80665;    % acceleration of gravity
    gradISA     = -0.0065; % Temperature gradient versus height K/m 
    Temp        = 20; %°C
    Pressure    = 1; %atm
    DegToRad    = pi/180;
    RadToDeg    = DegToRad^-1;
    msTokmh     = 60^2/1000;
    
	% b) Vehicle geometry and parameters (Ellipsoid) [m,Kg]
    a       = 0.9;       % x
    b       = 0.45;      % y
    c       = 0.45;      % z
    Vol     = 4*pi*a*b*c/3;
    Area    = Vol^(2/3);
    mBlimp  = Vol*(rhoAir-rhoHe);	% Mass
    mHe     = rhoHe*Vol;	% Mass due to Helium
    mTotal  = mBlimp+mHe;	% Total mass
    dx      = a/8;	% X axis distance from CV to propellers
    dy      = b/2;  % Y axis distance from CV to propellers
    dz      = c;    % Z axis distance from CV to propellers
    ax      = 0;    % X axis distance from center of gravity CG and center of volume CV
    ay      = 0;    % Z axis distance from center of gravity CG and center of volume CV
    az      = -(0.2*mBlimp)*b/mTotal;  % Z axis distance from center of gravity CG and center of volume CV
    displacedAirMass = Vol*rhoAir;

    % c) Masses and inertias
    Ix = mTotal*(b^2+c^2)/5; % Inertia of an ellipsoid along the x axis
    Iy = mTotal*(c^2+a^2)/5; % Inertia of an ellipsoid along the y axis
    Iz = mTotal*(a^2+b^2)/5; % Inertia of an ellipsoid along the z axis
    
    % c.1) Tuckerman for a prolate ellipsoid
    eTuckerman = sqrt(1-c^2/a^2);
    alphaTuckerman = (1-eTuckerman^2)*(log((1+eTuckerman)/(1-eTuckerman))-2*eTuckerman)/eTuckerman^3;
    betaTuckerman = (1-eTuckerman^2)*((eTuckerman/(1-eTuckerman^2))-0.5*log((1+eTuckerman)/(1-eTuckerman)))/eTuckerman^3;
    gammaTuckerman = betaTuckerman;
    
    K1 = Vol*(alphaTuckerman/(2-alphaTuckerman));
    K2 = Vol*(betaTuckerman/(2-betaTuckerman));
    K3 = Vol*(gammaTuckerman/(2-gammaTuckerman));
    K1_ = Vol*Ix*(((b^2-c^2)/(b^2+c^2))^2*((gammaTuckerman-betaTuckerman)/(2*((b^2-c^2)/(b^2+c^2))-(gammaTuckerman-betaTuckerman))));
    K2_ = Vol*Iy*(((c^2-a^2)/(c^2+a^2))^2*((alphaTuckerman-gammaTuckerman)/(2*((c^2-a^2)/(c^2+a^2))-(alphaTuckerman-gammaTuckerman))));
    K3_ = Vol*Iz*(((a^2-b^2)/(a^2+b^2))^2*((betaTuckerman-alphaTuckerman)/(2*((a^2-b^2)/(a^2+b^2))-(betaTuckerman-alphaTuckerman))));
    
    % c.2) Virtual masses and inertias
        % Tuckerman
            Xu = -K1*rhoAir;
            Yv = -K2*rhoAir;
            Zw = -K3*rhoAir;
            Lp = 0;
            Mq = -K2_*rhoAir;
            Nr = -K3_*rhoAir;
        % Others (Gomes)
            Mu = 0;
            Lv = 0;
            Nv = 0;
            Mw = 0;
            Yp = 0;
            Xq = 0;
            Zq = 0;
            Yr = 0;
        % Group
            mx = mTotal-Xu;
            my = mTotal-Yv;
            mz = mTotal-Zw;
            Jx = Ix-Lp;
            Jy = Iy-Mq;
            Jz = Iz-Nr;
            Jxz = 0;

    % d) M matrix
        M = [mx             0                   0               0               mTotal*az-Xq    0;...
             0              my                  0               -mTotal*az-Yp   0               mTotal*ax-Yr;...
             0              0                   mz              0               -mTotal*ax-Zq   0;...
             0              -mTotal*az-Lv       0               Ix-Lp           0               -Jxz;...
             mTotal*az-Mu   0                   -mTotal*ax-Mw   0               Iy-Mq           0;...
             0              mTotal*ax-Nv         0               -Jxz            0               Iz-Nr];
     
        invM = inv(M); % M matrix inverse
        
%% Step 2: Simulation configuration
fprintf('\n\nStep 2: Simulation configuration...');
    % a) Time definition
    fprintf('\n\ta) Time definition...');
    ti = 10.1;
    step = 0.1;
    tf = 560;
    tt = ti:step:tf;
    nData = length(tt);
    
    % b) Process configuration and declaration
    fprintf('\n\tb) Process declaration...');
    mInputs = 3;
    mStates = 6;
    tol     = 1e-6;
    
    UT= zeros(mInputs,nData);
    UT(2,:)= 0.5*ones(1,nData);
    U = zeros(mInputs,nData);
    X = zeros(mStates,nData);   % Zero initial conditions
    Y = zeros(mStates,nData);   % Zero initial conditions
    Vearth = zeros(3,nData);    % Zero initial conditions
    Xearth = zeros(3,nData);    % Zero initial conditions
    DCM = zeros(3);             % DCM
    lambda = zeros(3);          % Lambda's
    dX = zeros(mStates,nData);  % Zero initial conditions
    Vtot = zeros(1,nData);      % Rectiliniar speed
    wTot = zeros(1,nData);      % Angular speed
    G  = zeros(mStates,nData);  % Gravitational force
    A  = zeros(mStates,nData);  % Aerodynamic force
    Fd = zeros(mStates,nData);  % Dynamic forces
    P  = zeros(mStates,nData);  % Propulsion forces
    
    % c) Tools
    fprintf('\n\tc) Signal tools...');
    h = @(t) heaviside(t);
    r = @(t) heaviside(t)*t;    
    
    % d) Input generation
    fprintf('\n\td) Generating inputs...');
    Amp1  = 0.1;            % Motor maximum forces [N]
    Amp2  = 0.05;           % Motor thrust difference
    Amp3 = 90*DegToRad;     % Motor maximum angle rotation mu
    u1 = @(t) h(t)*-Amp1+r(t-10)*Amp1*0.0075-r(t-280)*Amp1*0.0075-(+r(t-290)*Amp1*0.0075-r(t-560)*Amp1*0.0075);
    u2 = @(t) h(t)*0.5+r(t-10)*Amp2*0.015-r(t-100)*Amp2*0.03+r(t-280)*Amp2*0.03-r(t-450)*Amp2*0.03;
    u3 = @(t) h(t)*-Amp3+r(t-10)*Amp3*0.0075-r(t-280)*Amp3*0.0075-h(t-290)*Amp3*2+r(t-290)*Amp3*0.0075-r(t-560)*Amp3*0.0075; 
    for i=1:nData
        UT(1,i) = (u1(tt(i)))+(randn-0.5)*Amp1*0.1;
        UT(2,i) = (u2(tt(i)))+(randn-0.5)*Amp2*0.1;
        UT(3,i) = (u3(tt(i)))+(randn-0.5)*Amp3*0.1;
        if UT(1,i)<-Amp1
            UT(1,i)=-Amp1;
        end
        if UT(1,i)>Amp1
            UT(1,i)=Amp1;
        end
        if UT(2,i)<0.5-Amp2
            UT(2,i)=0.5-Amp2;
        end
        if UT(2,i)>0.5+Amp2
            UT(2,i)=0.5+Amp2;
        end
        if UT(3,i)<-Amp3
            UT(3,i)=-Amp3;
        end
        if UT(3,i)>Amp3
            UT(3,i)=Amp3;
        end
    end


%% Step 3: Simulation 
fprintf('\n\nStep 3: Blimp simulation start...');
  for n=2:nData
    % a) Dynamics vector, Fd
    f1 = -mz*X(3,n-1)*X(5,n-1)+my*X(6,n-1)*X(2,n-1)+mTotal*(ax*(X(5,n-1)^2+X(6,n-1)^2)-az*X(6,n-1)*X(4,n-1));
    f2 = -mx*X(1,n-1)*X(6,n-1)+mz*X(4,n-1)*X(3,n-1)+mTotal*(-ax*X(4,n-1)*X(5,n-1)-az*X(6,n-1)*X(5,n-1));
    f3 = -my*X(2,n-1)*X(4,n-1)+mx*X(5,n-1)*X(1,n-1)+mTotal*(-ax*X(6,n-1)*X(4,n-1)+az*(X(5,n-1)^2+X(4,n-1)^2));
    f4 = -(Jz-Jy)*X(6,n-1)*X(5,n-1)+Jxz*X(4,n-1)*X(5,n-1)+mTotal*az*(X(1,n-1)*X(6,n-1)-X(4,n-1)*X(3,n-1));
    f5 = -(Jx-Jz)*X(4,n-1)*X(6,n-1)+Jxz*(X(6,n-1)^2-X(4,n-1)^2)+mTotal*(ax*(X(2,n-1)*X(4,n-1)-X(5,n-1)*X(1,n-1))-az*(X(3,n-1)*X(5,n-1)-X(6,n-1)*X(2,n-1)));
    f6 = -(Jy-Jx)*X(5,n-1)*X(4,n-1)-Jxz*X(5,n-1)*X(6,n-1)+mTotal*(-ax*(X(1,n-1)*X(6,n-1)-X(4,n-1)*X(3,n-1)));
    Fd(:,n-1) = [f1 f2 f3 f4 f5 f6]';
    
    % Input transforming   
    U(1,n-1)= UT(1,n-1)*UT(2,n-1); % Alpha*Tmax
    U(2,n-1)= UT(1,n-1)*(1-UT(2,n-1)); % (1-Alpha)*Tmax
    U(3,n-1)= UT(3,n-1);
    
    % b) Propulsion vector, P
    P1 = (U(1,n-1)+U(2,n-1))*cos(U(3,n-1));
    P2 = 0;
    P3 = -(U(1,n-1)+U(2,n-1))*sin(U(3,n-1));
    P4 = (U(2,n-1)-U(1,n-1))*sin(U(3,n-1))*dy;
    P5 = (U(1,n-1)+U(2,n-1))*(dz*cos(U(3,n-1))-dx*sin(U(3,n-1)));
    P6 = (U(2,n-1)-U(1,n-1))*cos(U(3,n-1))*dy; 
    P(:,n-1) =[P1 P2 P3 P4 P5 P6]';
    
    
    % c) Aerodynamic force vector, A   
    % c.1) Absolute speeds
      Vtot(n-1) = sqrt(X(1,n-1)^2+X(2,n-1)^2+X(3,n-1)^2);
      wTot(n-1) = sqrt(X(4,n-1)^2+X(5,n-1)^2+X(6,n-1)^2);
    
    % c.2) Aerodynamic coefficients
      CD = 0.9;
      CY = 0.9;
      CL = 0.9;
      Cl = 0.9;
      Cm = 0.9;
      Cn = 0.9; 
        
      coefA1 = 0.5*rhoAir*Vtot(n-1)^2*Area;
      coefA2 = 0.5*rhoAir*Vtot(n-1)^2*Vol;
      coefB1 = 0.5*rhoAir*X(1,n-1)^2*sign(X(1,n-1))*Area;
      coefB2 = 0.5*rhoAir*X(2,n-1)^2*sign(X(2,n-1))*Area;
      coefB3 = 0.5*rhoAir*X(3,n-1)^2*sign(X(3,n-1))*Area;
      coefB4 = 0.5*rhoAir*X(4,n-1)^2*sign(X(4,n-1))*Vol;
      coefB5 = 0.5*rhoAir*X(5,n-1)^2*sign(X(5,n-1))*Vol;
      coefB6 = 0.5*rhoAir*X(6,n-1)^2*sign(X(6,n-1))*Vol;
    
    % c.3) Aerodynamic forces A
      A1 = -CD*coefB1;
      A2 = -CY*coefB2;
      A3 = -CL*coefB3;
      A4 = -Cl*coefB4;
      A5 = -Cm*coefB5;
      A6 = -Cn*coefB6;
      A(:,n-1) = [A1 A2 A3 A4 A5 A6]';
    
    % d) Gravitational force vector, G 
    % d.1) Calculate the direction cosines matrix (DCM)
      lambda(1,1) = cos(Y(5,n-1))*cos(Y(6,n-1));
      lambda(1,2) = cos(Y(5,n-1))*sin(Y(6,n-1));
      lambda(1,3) = sin(Y(5,n-1));
      lambda(2,1) = (-cos(Y(4,n-1))*sin(Y(6,n-1))+sin(Y(4,n-1))*sin(Y(5,n-1))*cos(Y(6,n-1)));
      lambda(2,2) = (cos(Y(4,n-1))*cos(Y(6,n-1))+sin(Y(4,n-1))*sin(Y(5,n-1))*sin(Y(6,n-1)));
      lambda(2,3) = sin(Y(4,n-1))*cos(Y(5,n-1));
      lambda(3,1) = (sin(Y(4,n-1))*sin(Y(6,n-1))+cos(Y(4,n-1))*sin(Y(5,n-1))*cos(Y(6,n-1)));
      lambda(3,2) = (-sin(Y(4,n-1))*cos(Y(6,n-1))+cos(Y(4,n-1))*sin(Y(5,n-1))*sin(Y(6,n-1)));
      lambda(3,3) = cos(Y(4,n-1))*cos(Y(5,n-1));
      DCM = lambda;        
    
    % d.2) Calculate gravitational forces & moments. 
      B = rhoAir*gacc*Vol; % constant
      W = mTotal*gacc;     % constant
    
	% d.3) Using DCM elements convert these to body axes to obtain gravity vector G
      G1 = lambda(3,1)*(W-B);
      G2 = lambda(3,2)*(W-B);
      G3 = lambda(3,3)*(W-B);
      G4 = -lambda(3,2)*az*W;
      G5 = (lambda(3,1)*az-lambda(3,3)*ax)*W;
      G6 = lambda(3,2)*ax*W;
      G(:,n-1) = [G1 G2 G3 G4 G5 G6]';
        
    % e) Calculate linear and angular accelerations du, dv, dw, dp, dq, dr from [acc] =M-1[[P]+[Fd]+[A]+[G]]
      dX(:,n-1) = invM*(P(:,n-1)+Fd(:,n-1)+A(:,n-1)+G(:,n-1));
    
    % f) Integrate accelerations and obtain linear and angular body axes velocities u, v, w, p, q, r
      for i=1:mStates
          X(i,n) = sum(dX(i,1:n-1))*step+step*(dX(i,n-1)-dX(i,1))/2;
          if abs(X(i,n))<tol
              %X(i,n) = tol;
          end
      end
    
    % g) Transform linear velocities to earth axes to obtain Vnorth, Veast, Vup
      Vearth(:,n) = DCM*X(1:3,n);
      for i=1:3
          Xearth(i,n) = sum(Vearth(i,1:n-1))*step+step*(Vearth(i,n-1)-Vearth(i,1))/2;
      end
    
    % h) Calculate vehicle position in terms of displacements in the north, east and vertical directions
      for i=1:mStates
          Y(i,n) = sum(X(i,1:n-1))*step+step*(X(i,n-1)-X(i,1))/2; % (XYZ) aligned with (NEU)
      end
    
    % end
    %fprintf('\nn = %d',n);
  end  
    
%% Step 4: Simulation results
fprintf('\n\nStep 4: Plotting simulation results...');
lineWidth = 1.5;

for Step4=1:1
    
% a) Figure(1): Inputs
    figure(1)
    subplot(5,2,1);
    cla
    grid minor
    hold on;
    plot(tt,U(1,:),'-','LineWidth',lineWidth);
    legend('u1');
    title('Starboard motor force Tds'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,2,3);
    cla
    grid minor
    hold on;
    plot(tt,U(2,:),'-','LineWidth',lineWidth);
    legend('u2');
    title('Port motor force Tsp'); %xlabel('Time (s)'); 
    ylabel('(N)');
    
    subplot(5,2,2);
    cla
    grid minor
    hold on;
    plot(tt,UT(1,:),'-','LineWidth',lineWidth);
    legend('uT1');
    title('Motor max thrust Tmax'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,2,4);
    cla
    grid minor
    hold on;
    plot(tt,UT(2,:),'-','LineWidth',lineWidth);
    legend('uT2');
    title('Motor thrust difference ALPHA'); %xlabel('Time (s)'); 
    ylabel('[]');

    subplot(5,2,[5 6]);
    cla
    grid minor
    hold on;
    plot(tt,UT(3,:).*RadToDeg,'-','LineWidth',lineWidth);
    legend('u3');
    title('Propeller angle \mu'); %xlabel('Time (s)');  
    ylabel('(°)');

    subplot(5,2,[7 8]);
    cla
    grid minor
    hold on;
    plot(tt,Vtot,'-','LineWidth',lineWidth);
    %legend('');
    title('Net speed'); xlabel('Time (s)');  ylabel('(m/s)');

    subplot(5,2,[9 10]);
    cla
    grid minor
    hold on;
    plot(tt,wTot.*RadToDeg,'-','LineWidth',lineWidth);
    %legend('');
    title('Net angular speed'); xlabel('Time (s)');  ylabel('(°/s)');

% b) Figure(2): Velocities
    figure(2)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(tt,X(1,:),'-','LineWidth',lineWidth);
    title('u = X axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(tt,X(2,:),'-','LineWidth',lineWidth);
    title('v = Y axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(tt,X(3,:),'-','LineWidth',lineWidth);
    title('w = Z axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(tt,X(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('p = X axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(tt,X(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('q = Y axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(tt,X(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('r = Z axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
% c) Figure(3): Position and orientation
    figure(3)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(tt,Y(1,:),'-','LineWidth',lineWidth);
    title('x = X axis Position'); xlabel('Time (s)'); ylabel('Position (m)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(tt,Y(2,:),'-','LineWidth',lineWidth);
    title('y = Y axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(tt,Y(3,:),'-','LineWidth',lineWidth);
    title('z = Z axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(tt,Y(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\phi = X axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(tt,Y(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\theta = Y axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(tt,Y(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\psi = Z axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');

% d) Figure(4): Forces
    figure(4)
    subplot(2,2,1);
    cla
    grid minor
    hold on;
    plot(tt,Fd(1,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(2,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(3,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(4,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(5,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(6,:),'-','LineWidth',lineWidth);
    title('Fd(X(t))'); xlabel('Time (s)'); ylabel('Dynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,2);
    cla
    grid minor
    hold on;
    plot(tt,G(1,:),'-','LineWidth',lineWidth);
    plot(tt,G(2,:),'-','LineWidth',lineWidth);
    plot(tt,G(3,:),'-','LineWidth',lineWidth);
    plot(tt,G(4,:),'-','LineWidth',lineWidth);
    plot(tt,G(5,:),'-','LineWidth',lineWidth);
    plot(tt,G(6,:),'-','LineWidth',lineWidth);
    title('G(X(t))'); xlabel('Time (s)'); ylabel('Gravity force (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,3);
    cla
    grid minor
    hold on;
    plot(tt,A(1,:),'-','LineWidth',lineWidth);
    plot(tt,A(2,:),'-','LineWidth',lineWidth);
    plot(tt,A(3,:),'-','LineWidth',lineWidth);
    plot(tt,A(4,:),'-','LineWidth',lineWidth);
    plot(tt,A(5,:),'-','LineWidth',lineWidth);
    plot(tt,A(6,:),'-','LineWidth',lineWidth);
    title('A(X(t))'); xlabel('Time (s)'); ylabel('Aerodynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,4);
    cla
    grid minor
    hold on;
    plot(tt,P(1,:),'-','LineWidth',lineWidth);
    plot(tt,P(2,:),'-','LineWidth',lineWidth);
    plot(tt,P(3,:),'-','LineWidth',lineWidth);
    plot(tt,P(4,:),'-','LineWidth',lineWidth);
    plot(tt,P(5,:),'-','LineWidth',lineWidth);
    plot(tt,P(6,:),'-','LineWidth',lineWidth);
    title('P(u(t))'); xlabel('Time (s)'); ylabel('Propulsion forces (N)');
    legend('x','y','z','roll','pitch','yaw');  
    
end
    
%% Step 5: Scaling results to fuzzy sets
fprintf('\n\nStep 5: Scaling results to fuzzy sets...');
Mapping = @(Input,minIn,maxIn,minOut,maxOut) ((Input-minIn)./(maxIn-minIn)).*(maxOut-minOut)+minOut;

% Parámetros
minOut = 10.0;  % Porcentaje
maxOut = 90.0;  % Porcentaje
minUT1 = -0.1;  % Fuerza (N)
maxUT1 =  0.1;  % Fuerza (N)
minUT2 =  0.45; % Adimensional
maxUT2 =  0.55; % Adimensional
minUT3 = -pi/2; % Ángulo (rad)
maxUT3 =  pi/2; % Ángulo (rad)
minX1 = -1.0;   % Velocidad (m/s)
maxX1 =  1.0;   % Velocidad (m/s)
minX2 = -1.0;   % Velocidad (m/s)
maxX2 =  1.0;   % Velocidad (m/s)
minX3 = -1.0;   % Velocidad (m/s)
maxX3 =  1.0;   % Velocidad (m/s)
minX4 = -pi;    % Velocidad angular (rad/s)
maxX4 =  pi;    % Velocidad angular (rad/s)
minX5 = -pi;    % Velocidad angular (rad/s)
maxX5 =  pi;    % Velocidad angular (rad/s)
minX6 = -pi;    % Velocidad angular (rad/s)
maxX6 =  pi;    % Velocidad angular (rad/s)
minY1 = -100;   % Posición (m)
maxY1 =  100;   % Posición (m)
minY2 = -100;   % Posición (m)
maxY2 =  100;   % Posición (m)
minY3 = -100;   % Posición (m)
maxY3 =  100;   % Posición (m)
minY4 = -pi;    % Orientación (rad)
maxY4 =  pi;    % Orientación (rad)
minY5 = -pi;    % Orientación (rad)
maxY5 =  pi;    % Orientación (rad)
minY6 = -pi;    % Orientación (rad)
maxY6 =  pi;    % Orientación (rad)

% Conversión
uuu = zeros(mInputs,nData);
xxx = zeros(mStates,nData);
yyy = zeros(mStates,nData);

uuu(1,:) = Mapping(UT(1,:),minUT1,maxUT1,minOut,maxOut);
uuu(2,:) = Mapping(UT(2,:),minUT2,maxUT2,minOut,maxOut);
uuu(3,:) = Mapping(UT(3,:),minUT3,maxUT3,minOut,maxOut);

xxx(1,:) = Mapping(X(1,:),minX1,maxX1,minOut,maxOut);
xxx(2,:) = Mapping(X(2,:),minX2,maxX2,minOut,maxOut);
xxx(3,:) = Mapping(X(3,:),minX3,maxX3,minOut,maxOut);
xxx(4,:) = Mapping(X(4,:),minX4,maxX4,minOut,maxOut);
xxx(5,:) = Mapping(X(5,:),minX5,maxX5,minOut,maxOut);
xxx(6,:) = Mapping(X(6,:),minX6,maxX6,minOut,maxOut);

yyy(1,:) = Mapping(Y(1,:),minY1,maxY1,minOut,maxOut);
yyy(2,:) = Mapping(Y(2,:),minY2,maxY2,minOut,maxOut);
yyy(3,:) = Mapping(Y(3,:),minY3,maxY3,minOut,maxOut);
yyy(4,:) = Mapping(Y(4,:),minY4,maxY4,minOut,maxOut);
yyy(5,:) = Mapping(Y(5,:),minY5,maxY5,minOut,maxOut);
yyy(6,:) = Mapping(Y(6,:),minY6,maxY6,minOut,maxOut);
    
%% Step 6: Training configuration and space generation
fprintf('\n\nStep 6: ANFIS training configuration adn space generation...');
  
% a) Training configuration
  fprintf('\n\ta) Training configuration...');
  xx = 0:0.01:100;              % For plotting MFs
  nTrainingData = nData;        
  maxEpochs = 1000;            	% The more the merrier
  nInputs = 4;                  % Number of inputs, needs to be configured an ANFIS
  nFuzzy = 5;                   % Number of MFs per input
  nOutputs = 1;                 % Not changeable
  nRules = nFuzzy^nInputs;      % Number of rules
  K = 0.02;                     % Initial K
  maxK = 0.25;                 	% Maximum K, doesn't pass this number
  growthRate = 0.1;             % Growth of K 
  B = 1;                        % Backward values
  aIno = 15;                   	% Initial a parameter of premise MFs
  aOuto = 0.01;                 % Initial a parameter of consequent MFs
  tol = 1e-6;                 	% Tolerance for divisions by 0
  
% b) I/O
  fprintf('\n\tb) Selecting Inputs and outputs...');
  OUTPUT = xxx(1,:);      % Target function
  INPUT1 = uuu(1,:);      % First input        
  INPUT2 = uuu(2,:);      % Second input
  INPUT3 = uuu(3,:);      % Third input
  INPUT4 = xxx(1,:);      % Fourth input
  INPUT5 = xxx(2,:);      % Fifth input
  INPUT6 = xxx(3,:);      % Sixth input
  INPUT7 = yyy(3,:);      % Sixth input
  INPUT8 = uuu(1,:);      % Sixth input
  INPUT9 = uuu(1,:);      % Sixth input
  
% c) Tools
  fprintf('\n\tc) Generating tools...');
  gauss = @(x,a,c) exp(-(x-c).^2./(a.^2));      % Gauss MF for premise part
  sigmoid = @(x,a,c) 1./(1+exp(-(x-c)./a));     % Sigmoid MF for consequent par

  invSigmoid = @(x,a,c) c-a.*log(1./x-1);       % Inverse of sigmoid function

  dGauss_da = @(x,a,c) (2.*exp(-(-c+x).*(-c+x)./(a.*a)).*(-c+x).*(-c+x))./(a.*a.*a);    % Partial w/r to a
  dGauss_dc = @(x,a,c) (2.*exp(-(-c+x).*(-c+x)./(a.*a)).*(-c+x))./(a.*a);               % Partial w/r to c

  dinvSigmoid_da = @(x,a,c) -log(1./x-1);   % Partial of invSigmoid w/r to a
  dinvSigmoid_dc = @(x,a,c) 1;              % Partial of invSigmoid w/r to c

% d) Generating workspace
 fprintf('\n\td) Generating workspace...');
  % d.1) Initial fuzzy parameters
    fprintf('\n\t\td.1) Initial fuzzy parameters...');
    aIne  = zeros(nFuzzy,nInputs);
    cIne  = zeros(nFuzzy,nInputs);
    aOute = zeros(nRules,nOutputs);
    cOute = zeros(nRules,nOutputs);
    aInfinal  = aIne;
    cInfinal  = cIne;
    aOutfinal = aOute;
    cOutfinal = cOute;
    for i=1:nInputs
        aIne(:,i) = aIno.*ones(nFuzzy,1);   % Initial a for premise MFs
        cIne(:,i) = (0:(100/(nFuzzy-1)):100); % Initial c for premise MFs
    end
    for i=1:nOutputs
        aOute(:,i) = aOuto*ones(nRules,1);   % Initial a for consequent MFs
        cOute(:,i) = (0:(100/(nRules-1)):100); % Initial c for consequent MFs
    end

  % d.2) Training workspace
    fprintf('\n\t\td.2) Training workspace...');
    APE         = zeros(maxEpochs,1);
    APEmin      = 100000;
    epochFlag   = 1;
    etaIna      = zeros(nFuzzy,nInputs);
    etaInc      = zeros(nFuzzy,nInputs);
    etaOuta     = zeros(nRules,1);
    etaOutc     = zeros(nRules,1);
    X           = zeros(nInputs,1);
    O5          = zeros(nTrainingData,1);
    En          = zeros(nTrainingData,1);
    muIn        = zeros(nFuzzy,nInputs);
    w           = zeros(nRules,1);
    wn          = zeros(nRules,1);
    fi          = zeros(nRules,1);
    fi_wn       = zeros(nRules,1);

    dJn_dO5     = zeros(nTrainingData,1);
    dJn_dO2     = zeros(nRules,1);
    dO5_dfi     = zeros(nRules,1);
    dO5_dO2     = zeros(1,nRules);
    dO2_dO1     = zeros(nRules,nFuzzy*nInputs);

    dfi_da      = zeros(nRules,1);
    dfi_dc      = zeros(nRules,1);
    dmu_daIn    = zeros(nFuzzy,nInputs);
    dmu_dcIn    = zeros(nFuzzy,nInputs);
    dJn_daOut   = zeros(nRules,1);
    dJn_dcOut   = zeros(nRules,1);
    dJp_daOut   = zeros(nRules,1);
    dJp_dcOut   = zeros(nRules,1);
    dJn_dmu     = zeros(nFuzzy,nInputs);
    dJn_daIn    = zeros(nFuzzy,nInputs);
    dJn_dcIn    = zeros(nFuzzy,nInputs);
    dJp_daIn    = zeros(nFuzzy,nInputs);
    dJp_dcIn    = zeros(nFuzzy,nInputs);
    
  % d.3) Index table
    fprintf('\n\t\td.3) Index table...');
    indexTable  = zeros(nRules,nInputs); 
    for k = 1:nInputs
        l=1;
        for j=1:nFuzzy^(k-1):nRules
            for i =1:nFuzzy^(k-1)
                indexTable(j+i-1,nInputs-k+1)=l;
            end
            l=l+1;
            if l>nFuzzy
                l=1;
            end
        end
    end 

% e) Plotting part II   
  fprintf('\n\te) Plotting part II...');
  figure(5)
  for i=1:nInputs
     subplot(nOutputs+1,nInputs,i); 
     cla
     hold on
     for j=1:nFuzzy
        plot(xx,gauss(xx,aIne(j,i),cIne(j,i)),'DisplayName',[num2str(char(64+i)) num2str(j)],'LineWidth',0.5);
     end
     %legend;
     title(['Initial premise fuzzy set ' num2str(char(64+i))]); ylabel(['\mu(x_',num2str(i),')']); xlabel(['x_',num2str(i)]);
     grid minor;
  end
  
  for i=1:nOutputs
    subplot(nOutputs+1,nInputs,[(i*nInputs+1) (i*nInputs+nInputs)]);
    cla
    hold on;
    for j=1:nRules
        plot(xx,sigmoid(xx,aOute(j,i),cOute(j,i)),'LineWidth',0.1);
    end
    %legend;
    title(['Initial consequent fuzzy set Z_' num2str(i)]); ylabel(['\mu(z_',num2str(i),')']); xlabel(['z_',num2str(i)]);
    grid minor
  end 
    
%% Step 7: ANFIS offline training
fprintf('\n\nStep 7: ANFIS offline training...');

for g=1:maxEpochs
  dJp_daOut	= zeros(nRules,1);
  dJp_dcOut	= zeros(nRules,1);
  dJp_daIn  = zeros(nFuzzy,nInputs);
  dJp_dcIn	= zeros(nFuzzy,nInputs); 
	for i=1+B:nTrainingData
	  %% a) ANFIS;
	  % Prelayer: Selecting input variables
	    X(1)  = INPUT1(i-B);
        X(2)  = INPUT2(i-B);
        X(3)  = INPUT3(i-B);
        X(4)  = INPUT4(i-B);
        %X(5)  = INPUT5(i-B);
        %X(6)  = INPUT6(i-B);
        
	  % Layer 1: Input fuzzyfication  
        for k=1:nInputs
            for j=1:nFuzzy
                muIn(j,k) = gauss(X(k),aIne(j,k),cIne(j,k));
                if muIn(j,k)<tol
                   muIn(j,k) = tol;
                end
                if muIn(j,k)>(1.0-tol)
                   muIn(j,k) = 1.0-tol;
                end
            end            
        end

	  % Layer 2: Calculating weigths
        for j=1:nRules
                w(j) = 1.0;
            for k=1:nInputs
                %w(j) = min(w(j),muIn(indexTable(j,k),k)); % Not recommendable
                w(j) = w(j)*muIn(indexTable(j,k),k);
            end
        end
        
	  % Layer 3: Normalizing
        sumW = sum(w);
        if abs(sumW)<tol
            sumW=tol;
        end
        for j=1:nRules
           wn(j) = w(j)/sumW;
        end
       
	  % Layer 4: Calculating wn_i*f_i     
        for j=1:nRules
           fi(j)   = invSigmoid(w(j),aOute(j,1),cOute(j,1));
           fi_wn(j) = fi(j)*wn(j); 
        end

	  % Layer 5: Calculating output
        f = sum(fi_wn);
        O5(i) = f;

	  %% b) ANFIS error measure
        En(i) = OUTPUT(i)-O5(i);      % Measured error for the i-th data pair
	  % dJn/dO5
	    dJn_dO5(i) = -2*En(i);
         
      %% c) Gradient descent for consequent parameters
        for j=1:nRules
        %dO5_dfi
          dO5_dfi(j) = wn(j);
        %dfi_dalpha
          dfi_da(j) = dinvSigmoid_da(w(j),aOute(j,1),cOute(j,1));
          dfi_dc(j) = 1; %dinvSigmoid_dc(w(j),aOute(j,1),cOute(j,1));
        % dJn_dalpha
          dJn_daOut(j) = dJn_dO5(i)*dO5_dfi(j)*dfi_da(j);
          dJn_dcOut(j) = dJn_dO5(i)*dO5_dfi(j)*dfi_dc(j);
        % Sum
          dJp_daOut(j) = dJp_daOut(j)+dJn_daOut(j);
          dJp_dcOut(j) = dJp_dcOut(j)+dJn_dcOut(j);
        end     

      %% d) Gradient descent for premise parameters
	  % dO5/dO2       
        for j=1:nRules
          dO5_dO2(j) = (fi(j)-O5(i))/sumW; 
        end
	  % dO2_dO1 matrix
        for e=1:nInputs
            for k=1:nFuzzy
                for j=1:nRules
                    if(k==indexTable(j,e))
                        dO2_dO1(j,(e-1)*nFuzzy+k)=1.0;
                        for p=1:nInputs
                            if(e~=p)
                                %dO2_dO1(j,(e-1)*nFuzzy+k) = min(dO2_dO1(j,(e-1)*nFuzzy+k),muIn(indexTable(j,p),p)); % Not recommendable
                                dO2_dO1(j,(e-1)*nFuzzy+k) = dO2_dO1(j,(e-1)*nFuzzy+k)*muIn(indexTable(j,p),p);
                            end
                        end
                    else
                        dO2_dO1(j,(e-1)*nFuzzy+k)=0.0;
                    end
                end
            end
        end     
	  % dJn_dO2
      for j=1:nRules
        dJn_dO2(j) = dJn_dO5(i)*dO5_dO2(j);  
      end
        
	  % Chain rule
        for k=1:nInputs
          for j=1:nFuzzy
          % dJn_dO1
            dJn_dmu(j,k) = 0.0;
            for p=1:nRules
                dJn_dmu(j,k) = dJn_dmu(j,k)+ dJn_dO2(p)*dO2_dO1(p,j+(k-1)*nFuzzy);
            end
          % dO1_dalpha
            dmu_daIn(j,k)= dGauss_da(X(k),aIne(j,k),cIne(j,k));
            dmu_dcIn(j,k)= dGauss_dc(X(k),aIne(j,k),cIne(j,k));
          % dJn_dalpha
            dJn_daIn(j,k) = (dJn_dmu(j,k)).*dmu_daIn(j,k);
            dJn_dcIn(j,k) = (dJn_dmu(j,k)).*dmu_dcIn(j,k);
          % Sum
            dJp_daIn(j,k) = dJp_daIn(j,k)+dJn_daIn(j,k);
            dJp_dcIn(j,k) = dJp_dcIn(j,k)+dJn_dcIn(j,k);
          end
        end
    end
    
    %% e) Epoch summary
      APE(g,1) = 0;
      for i=1+B:nTrainingData
      if abs(OUTPUT(i))<1
          OUTPUT(i) = 1;
      end
          APE(g,1) = APE(g,1)+abs(En(i))/abs(OUTPUT(i));
      end
      APE(g,1) = APE(g,1)*100/(nTrainingData-B);    
      if APE(g,1)<=APEmin
        APEmin = APE(g,1);
        aInfinal = aIne;
        cInfinal = cIne;
        aOutfinal = aOute;
        cOutfinal = cOute;
        epochFlag = g;
      end
    
    %% f) New step size    
    if g>4
        if APE(g,1)<APE(g-1,1)
            if APE(g-1,1)<APE(g-2,1)
                if APE(g-2,1)<APE(g-3,1)
                    if APE(g-3,1)<APE(g-4,1)
                        K = K*(1.0+growthRate);
                    end
                end
            end
        else
            if APE(g-1,1)<APE(g-2,1)
                if APE(g-2,1)>APE(g-3,1)
                    if APE(g-3,1)<APE(g-4,1)
                    K=K*(1.0-growthRate);          
                    end
                end
            end
        end
    end
    if K>maxK
        K = maxK;
    end
    
    %% g) New consequent parameters
	for j=1:nRules
	  aOute(j,1) = aOute(j,1)-K*sign(dJp_daOut(j,1));
	  cOute(j,1) = cOute(j,1)-K*sign(dJp_dcOut(j,1));
      if abs(aOute(j,1))<tol
          aOute(j,1)=tol;
      end
      if abs(cOute(j,1))<tol
          cOute(j,1)=tol;
      end
    end
        
  %% h) New premise parameters
%     for k=1:nInputs
%       for j=1:nFuzzy
%         aIne(j,k) = aIne(j,k)-K*sign(dJp_daIn(j,k));
%         cIne(j,k) = cIne(j,k)-K*sign(dJp_dcIn(j,k));
%         if abs(aIne(j,k))<tol
%           aIne(j,k)=tol;
%         end
%         if abs(cIne(j,k))<tol
%             cIne(j,k)=tol;
%         end
%       end
%     end
    
  % i) Print epoch summary
   fprintf('\n\tEpoch: %04d.\t APE: %.10f.\t K: %.10f.\t APEmin: %d,\t %.7f',g,APE(g),K, epochFlag,APEmin);
end

%% Step 8: Plotting training results
fprintf('\n\nStep 8: Plotting training results...');
  
  % a) APE evolution
  fprintf('\n\ta) APE evolution...');
  figure(6)
  hold on;
  subplot(2,1,1); 
  cla;
  hold on
  plot(1:length(APE(:,1)), APEmin.*ones(length(APE(:,1)),1), '--k', 'LineWidth', 0.1);
  plot(1:length(APE(:,1)), APE(:,1), '-', 'LineWidth', 1);
  ylim([0 max(APE)]);
  title('Average percentage error (APE) evolution'); xlabel('Epochs'); ylabel('APE (%)');
  legend('min','APE')
  grid minor;
  
  % b) Target and actual
  fprintf('\n\tb) Target and actual...');
  subplot(2,1,2);
  cla;
  hold on;
  legend('off');
  plot(tt, O5, '-b', 'LineWidth', 1.2); 
  plot(tt, OUTPUT, '-g', 'LineWidth', 1.2);   
  legend('Calculated','Target');
  title('Anfis training results'); xlabel('Time (s)'); ylabel('(%)');
  grid minor;
  
  % c) Best fuzzy sets
  fprintf('\n\tc) Best fuzzy sets...');
  figure(7)
  for i=1:nInputs
     subplot(nOutputs+1,nInputs,i); 
     hold on
     for j=1:nFuzzy
        plot(xx,gauss(xx,aInfinal(j,i),cInfinal(j,i)),'DisplayName',[num2str(char(64+i)) num2str(j)],'LineWidth',0.5);
     end
     %legend;
     title(['Best premise fuzzy set ' num2str(char(64+i))]); ylabel(['\mu(x_',num2str(i),')']); xlabel(['x_',num2str(i)]);
     grid minor;
  end
  
  for i=1:nOutputs
    subplot(nOutputs+1,nInputs,[(i*nInputs+1) (i*nInputs+nInputs)]);
    hold on;
    for j=1:nRules
        plot(xx,sigmoid(xx,aOutfinal(j,i),cOutfinal(j,i)),'LineWidth',0.1);
    end
    %legend;
    title(['Best consequent fuzzy set Z_' num2str(i)]); ylabel(['\mu(z_',num2str(i),')']); xlabel(['z_',num2str(i)]);
    grid minor
  end 

  % d) Training summary
  fprintf('\n\td) Training summary...');
  figure(8);
  cla;
  hold on;
  plot(tt,INPUT1, '-', 'LineWidth', 0.1,'Color',[0.1,0.1,0.1]); 
  plot(tt,INPUT2, '-', 'LineWidth', 0.1,'Color',[0.2,0.2,0.2]);  
  plot(tt,INPUT3, '-', 'LineWidth', 0.1,'Color',[0.3,0.3,0.3]); 
  plot(tt,[0,INPUT4(1:end-1)], '-', 'LineWidth', 0.1,'Color',[0.4,0.4,0.4]);
  plot(tt,O5, '-b', 'LineWidth', 1.1); 
  plot(tt,OUTPUT, '-g', 'LineWidth', 1.1);
  title('ANFIS training summary'); xlabel('Time (s)'); ylabel('(%)');
  grid minor; 
  legend('ANFIS input1','ANFIS input2','ANFIS input3','ANFIS input4','ANFIS output','Target');
  
  
  % e) Saving variables
  save('blimp_ANFIS_x0.mat','aInfinal','cInfinal','aOutfinal','cOutfinal');  
  save('blimp_ANFIS_settings_x0.mat'...
        ,'nInputs'      ,'nOutputs'     ,'mStates'      ,'mInputs'...
        ,'gauss'        ,'sigmoid'      ,'invSigmoid'   ,'h'...
        ,'tol'          ,'minOut'       ,'maxOut'       ,'r'       ...
        ,'indexTable'   ,'nFuzzy'       ,'nRules', 'Mapping'       ...
        ,'minUT1' ,'maxUT1' ,'minUT2' ,'maxUT2' ,'minUT3' ,'maxUT3'...
        ,'minX1'  ,'maxX1'  ,'minX2'  ,'maxX2'  ,'minX3'  ,'maxX3'...
        ,'minX4'  ,'maxX4'  ,'minX5'  ,'maxX5'  ,'minX6'  ,'maxX6'...
        ,'minY1'  ,'maxY1'  ,'minY2'  ,'maxY2'  ,'minY3'  ,'maxY3'...
        ,'minY4'  ,'maxY4'  ,'minY5'  ,'maxY5'  ,'minY6'  ,'maxY6'...
        ,'rhoAir'  ,'gacc'  ,'DegToRad'  ,'RadToDeg'  ,'msTokmh'...
        ,'Ix' ,'Iy' ,'Iz' ,'Xu' ,'Yv' ,'Zw' ,'Lp' ,'Mq' ,'Nr'...
        ,'Mu' ,'Lv' ,'Mw' ,'Yp' ,'Xq' ,'Zq' ,'Yr' ,'Jxz','mTotal'...
        ,'mx' ,'my' ,'mz' ,'Jx' ,'Jy' ,'Jz' ,'M' ,'invM'...
        ,'dx' ,'dy' ,'dz' ,'ax' ,'ay' ,'az', 'Area', 'Vol');

%% Step 8 and a half: Loading ANFIS
fprintf('\n\nStep 8 and a half: Loading ANFIS...');
  clear all;
  load('blimp_ANFIS_x0.mat');
  load('blimp_ANFIS_settings_x0.mat');

% Step 9: Control configuration
fprintf('\n\nStep 9: Control configuration...');

% a) Time definition
  fprintf('\n\ta) Time definition...');
  tii   = 0.1;
  step  = 0.1;
  tff   = 600;
  ttt   = tii:step:tff;
  nctrlData = length(ttt);
  
% b) Configuration
  fprintf('\n\tb) Configuration...');
  % b.1) DMC
  fprintf('\n\t\tb.1) Configuring DMC options...');
  alpha = 0.9;
  Np = 4;
  Nc = 3;
  Nu = 1;
  Q = 1*eye(Np);            % Error weight
  R = 1*eye(Nc*mInputs);   % Control action weight
  
  R(1:Nc,1:Nc)                  = 5*R(1:Nc,1:Nc);   % 1st Control action weight
  R(Nc+1:2*Nc,Nc+1:2*Nc)        = 1000*R(Nc+1:2*Nc,Nc+1:2*Nc);   % 2nd Control action weight
  R(2*Nc+1:3*Nc,2*Nc+1:3*Nc)    = 5*R(2*Nc+1:3*Nc,2*Nc+1:3*Nc);   % 3rd Control action weight
  
  % b.2) ANFIS
  fprintf('\n\t\tb.2) Configuring DMC options...');

% c) Workspace
  fprintf('\n\tc) Creating workspace...');
  % c.1) DMC
  fprintf('\n\t\tc.1) DMC workspace...');
  gi    = zeros(Np,1);
  Gij   = zeros(Np,mInputs);
  Wref  = zeros(Np,1);
  Yfree = zeros(Np,1);
  Ystep = zeros(Np,mInputs);
  
  % c.2) Process 
  fprintf('\n\t\tc.2) Process workspace...');
  UU        = zeros(mInputs,nctrlData);     % Zero U initial conditions
  UTT       = zeros(mInputs,nctrlData);     % Zero U initial conditions
  UTT(2,:)  = 0.5*ones(1,nctrlData);     % Zero U initial conditions
  XX        = zeros(mStates,nctrlData);     % Zero X initial conditions
  YY        = zeros(mStates,nctrlData);     % Zero Y initial conditions
  Vearth    = zeros(3,nctrlData);           % Zero X initial conditions
  Xearth    = zeros(3,nctrlData);           % Zero X initial conditions
  DCM       = zeros(3);                     % DCM
  lambda    = zeros(3);                     % Lambda's
  dX        = zeros(mStates,nctrlData);     % Zero initial conditions
  Vtot      = zeros(1,nctrlData);           % Rectiliniar speed
  wTot      = zeros(1,nctrlData);           % Angular speed
  G         = zeros(mStates,nctrlData);     % Gravitational force
  A         = zeros(mStates,nctrlData);     % Aerodynamic force
  Fd        = zeros(mStates,nctrlData);     % Dynamic forces
  P         = zeros(mStates,nctrlData);     % Propulsion forces

  % c.3) ANFIS 
  fprintf('\n\t\tc.3) ANFIS workspace...');
  X     = 50*ones(nInputs,1);
  ukk   = 50*ones(mInputs,nctrlData);
  dukk  = zeros(mInputs*Nc);
  xkk   = 50*ones(mStates,nctrlData);
  ykk   = 50*ones(mStates,nctrlData);
  yee   = 50*ones(1,nctrlData);
  
% d) Reference definition
  fprintf('\n\td) Reference definition...');
  reference = @(t) h(t)*50-h(t-30)*10+h(t-250)*10+h(t-350)*10;
  ref = zeros(nctrlData,1);
  ref(:,1) = reference(ttt);
  for i=1:nctrlData
     ref(i) = ref(i)+(rand-0.5)*0.01; 
  end
  
% e) Reference filter
  fprintf('\n\te) Filtering reference...');
  wref = zeros(nctrlData,1);
  wref(1) = ref(1);
  for i=2:nctrlData
     wref(i) = alpha*wref(i-1)+(1-alpha)*ref(i); 
  end  
  
% f) Reference plot 
  fprintf('\n\tf) Plotting part VI...');
  figure(9)
  cla
  hold on
  plot(ttt, ref);
  plot(ttt, wref);
  grid minor
  title('Control reference'); xlabel('Time (s)'); ylabel('RPM');
  legend('ref','wref');
  
% Step 10: Control sequence
fprintf('\n\nStep 10: Control sequence start...');

  for i=3:nctrlData-Np
fprintf('\nfor i = %d/%d - ukk(1,i-1) = %.3f - ukk(2,i-1) = %.3f - ukk(3,i-1) = %.3f',i,nctrlData-Np,ukk(1,i-1),ukk(2,i-1),ukk(3,i-1));
      % a) OUTPUT reading
        % Blimp (yk)
        for Simulation=1:1
            % a) Dynamics vector, Fd
            f1 = -mz*XX(3,i-1)*XX(5,i-1)+my*XX(6,i-1)*XX(2,i-1)+mTotal*(ax*(XX(5,i-1)^2+XX(6,i-1)^2)-az*XX(6,i-1)*XX(4,i-1));
            f2 = -mx*XX(1,i-1)*XX(6,i-1)+mz*XX(4,i-1)*XX(3,i-1)+mTotal*(-ax*XX(4,i-1)*XX(5,i-1)-az*XX(6,i-1)*XX(5,i-1));
            f3 = -my*XX(2,i-1)*XX(4,i-1)+mx*XX(5,i-1)*XX(1,i-1)+mTotal*(-ax*XX(6,i-1)*XX(4,i-1)+az*(XX(5,i-1)^2+XX(4,i-1)^2));
            f4 = -(Jz-Jy)*XX(6,i-1)*XX(5,i-1)+Jxz*XX(4,i-1)*XX(5,i-1)+mTotal*az*(XX(1,i-1)*XX(6,i-1)-XX(4,i-1)*XX(3,i-1));
            f5 = -(Jx-Jz)*XX(4,i-1)*XX(6,i-1)+Jxz*(XX(6,i-1)^2-XX(4,i-1)^2)+mTotal*(ax*(XX(2,i-1)*XX(4,i-1)-XX(5,i-1)*XX(1,i-1))-az*(XX(3,i-1)*XX(5,i-1)-XX(6,i-1)*XX(2,i-1)));
            f6 = -(Jy-Jx)*XX(5,i-1)*XX(4,i-1)-Jxz*XX(5,i-1)*XX(6,i-1)+mTotal*(-ax*(XX(1,i-1)*XX(6,i-1)-XX(4,i-1)*XX(3,i-1)));
            Fd(:,i-1) = [f1 f2 f3 f4 f5 f6]';
            
            % Input transforming   
            UU(1,i-1)= UTT(1,i-1)*UTT(2,i-1); % Alpha*Tmax
            UU(2,i-1)= UTT(1,i-1)*(1-UTT(2,i-1)); % (1-Alpha)*Tmax
            UU(3,i-1)= UTT(3,i-1);
            
            % b) Propulsion vector, P
            P1 = (UU(1,i-1)+UU(2,i-1))*cos(UU(3,i-1));
            P2 = 0;
            P3 = -(UU(1,i-1)+UU(2,i-1))*sin(UU(3,i-1));
            P4 = (UU(2,i-1)-UU(1,i-1))*sin(UU(3,i-1))*dy;
            P5 = (UU(1,i-1)+UU(2,i-1))*(dz*cos(UU(3,i-1))-dx*sin(UU(3,i-1)));
            P6 = (UU(2,i-1)-UU(1,i-1))*cos(UU(3,i-1))*dy; 
            P(:,i-1) =[P1 P2 P3 P4 P5 P6]';


            % c) Aerodynamic force vector, A   
            % c.1) Absolute speeds
              Vtot(i-1) = sqrt(XX(1,i-1)^2+XX(2,i-1)^2+XX(3,i-1)^2);
              wTot(i-1) = sqrt(XX(4,i-1)^2+XX(5,i-1)^2+XX(6,i-1)^2);

            % c.2) Aerodynamic coefficients
              CD = 0.47;
              CY = 0.47;
              CL = 0.47;
              Cl = 0.47;
              Cm = 0.47;
              Cn = 0.47; 

              coefA1 = 0.5*rhoAir*Vtot(i-1)^2*Area;
              coefA2 = 0.5*rhoAir*Vtot(i-1)^2*Vol;
              coefB1 = 0.5*rhoAir*XX(1,i-1)^2*sign(XX(1,i-1))*Area;
              coefB2 = 0.5*rhoAir*XX(2,i-1)^2*sign(XX(2,i-1))*Area;
              coefB3 = 0.5*rhoAir*XX(3,i-1)^2*sign(XX(3,i-1))*Area;
              coefB4 = 0.5*rhoAir*XX(4,i-1)^2*sign(XX(4,i-1))*Vol;
              coefB5 = 0.5*rhoAir*XX(5,i-1)^2*sign(XX(5,i-1))*Vol;
              coefB6 = 0.5*rhoAir*XX(6,i-1)^2*sign(XX(6,i-1))*Vol;

            % c.3) Aerodynamic forces A
              A1 = -CD*coefB1;
              A2 = -CY*coefB2;
              A3 = -CL*coefB3;
              A4 = -Cl*coefB4;
              A5 = -Cm*coefB5;
              A6 = -Cn*coefB6;
              A(:,i-1) = [A1 A2 A3 A4 A5 A6]';

            % d) Gravitational force vector, G 
            % d.1) Calculate the direction cosines matrix (DCM)
              lambda(1,1) = cos(YY(5,i-1))*cos(YY(6,i-1));
              lambda(1,2) = cos(YY(5,i-1))*sin(YY(6,i-1));
              lambda(1,3) = sin(YY(5,i-1));
              lambda(2,1) = (-cos(YY(4,i-1))*sin(YY(6,i-1))+sin(YY(4,i-1))*sin(YY(5,i-1))*cos(YY(6,i-1)));
              lambda(2,2) = (cos(YY(4,i-1))*cos(YY(6,i-1))+sin(YY(4,i-1))*sin(YY(5,i-1))*sin(YY(6,i-1)));
              lambda(2,3) = sin(YY(4,i-1))*cos(YY(5,i-1));
              lambda(3,1) = (sin(YY(4,i-1))*sin(YY(6,i-1))+cos(YY(4,i-1))*sin(YY(5,i-1))*cos(YY(6,i-1)));
              lambda(3,2) = (-sin(YY(4,i-1))*cos(YY(6,i-1))+cos(YY(4,i-1))*sin(YY(5,i-1))*sin(YY(6,i-1)));
              lambda(3,3) = cos(YY(4,i-1))*cos(YY(5,i-1));
              DCM = lambda;        

            % d.2) Calculate gravitational forces & moments. 
              B = rhoAir*gacc*Vol; % constant
              W = mTotal*gacc;     % constant

            % d.3) Using DCM elements convert these to body axes to obtain gravity vector G
              G1 = lambda(3,1)*(W-B);
              G2 = lambda(3,2)*(W-B);
              G3 = lambda(3,3)*(W-B);
              G4 = -lambda(3,2)*az*W;
              G5 = (lambda(3,1)*az-lambda(3,3)*ax)*W;
              G6 = lambda(3,2)*ax*W;
              G(:,i-1) = [G1 G2 G3 G4 G5 G6]';

            % e) Calculate linear and angular accelerations du, dv, dw, dp, dq, dr from [acc] =M-1[[P]+[Fd]+[A]+[G]]
              dX(:,i-1) = invM*(P(:,i-1)+Fd(:,i-1)+A(:,i-1)+G(:,i-1));

            % f) Integrate accelerations and obtain linear and angular body axes velocities u, v, w, p, q, r
              for n=1:mStates
                  XX(n,i) = sum(dX(n,1:i-1))*step+step*(dX(n,i-1)-dX(n,1))/2;
              end

            % g) Transform linear velocities to earth axes to obtain Vnorth, Veast, Vup
              Vearth(:,i) = DCM*XX(1:3,i);
              for n=1:3
                  Xearth(n,i) = sum(Vearth(n,1:i-1))*step+step*(Vearth(n,i-1)-Vearth(n,1))/2;
              end

            % h) Calculate vehicle position in terms of displacements in the north, east and vertical directions
              for n=1:mStates
                  YY(n,i) = sum(XX(n,1:i-1))*step+step*(XX(n,i-1)-XX(n,1))/2; % (XYZ) aligned with (NEU)
              end
        end
        
        % Scaling to fuzzy
        xkk(1,i) = Mapping(XX(1,i),minX1,maxX1,minOut,maxOut);
        xkk(2,i) = Mapping(XX(2,i),minX2,maxX2,minOut,maxOut);
        xkk(3,i) = Mapping(XX(3,i),minX3,maxX3,minOut,maxOut);
        xkk(4,i) = Mapping(XX(4,i),minX4,maxX4,minOut,maxOut);
        xkk(5,i) = Mapping(XX(5,i),minX5,maxX5,minOut,maxOut);
        xkk(6,i) = Mapping(XX(6,i),minX6,maxX6,minOut,maxOut);

        ykk(1,i) = Mapping(YY(1,i),minY1,maxY1,minOut,maxOut);
        ykk(2,i) = Mapping(YY(2,i),minY2,maxY2,minOut,maxOut);
        ykk(3,i) = Mapping(YY(3,i),minY3,maxY3,minOut,maxOut);
        ykk(4,i) = Mapping(YY(4,i),minY4,maxY4,minOut,maxOut);
        ykk(5,i) = Mapping(YY(5,i),minY5,maxY5,minOut,maxOut);
        ykk(6,i) = Mapping(YY(6,i),minY6,maxY6,minOut,maxOut);
           
      % b) ANFIS OUTPUT reading
            % Prelayer: Input variables
            X(1) = ukk(1,i-1);
            X(2) = ukk(2,i-1);
            X(3) = ukk(3,i-1);
            X(4) = xkk(1,i-1);
            % Layer 1: Fuzzification 
            for k=1:nInputs
                for j=1:nFuzzy
                    muIn(j,k) = gauss(X(k),aInfinal(j,k),cInfinal(j,k));    % Better results
                    if muIn(j,k)<tol
                       muIn(j,k) = tol;
                    end
                    if muIn(j,k)>1.0-tol
                       muIn(j,k) = 1.0-tol;
                    end
                end            
            end
            %Layer 2: Calculation of w 
            for j=1:nRules
                    w(j) = 1.0;
                for k=1:nInputs
                    w(j) = w(j).*muIn(indexTable(j,k),k);
                end
            end
            % Layer 3: Normalization
            sumW = sum(w);
            if abs(sumW)<tol
                sumW=tol;
            end
            for j=1:nRules
               wn(j) = w(j)/sumW;
            end
            % Layer 4: Cálculo de wn_i*f_i 
            for j=1:nRules
               fi(j)   = invSigmoid(w(j),aOutfinal(j,1),cOutfinal(j,1)); % Better results
               fi_wn(j) = fi(j)*wn(j); 
            end
            % Layer 5: Final layer
             f = sum(fi_wn);
             yee(i) = f;

      % Wref
        for j=1:Np
            Wref(j) = wref(i);
        end
      
      % Error
        ekk = wref(i)-xkk(1,i);
        
      % c) Prediction
      for PREDICTION=1:1
        % ANFIS (Yfree)
        Yfree(1) = xkk(1,i-1);
       	for z=1:Np-1
        	% Prelayer: Input variables
            X(1) = ukk(1,i-1);
            X(2) = ukk(2,i-1);
            X(3) = ukk(3,i-1);
            X(4) = Yfree(z);
            % Layer 1: Fuzzification 
            for k=1:nInputs
                for j=1:nFuzzy
                    muIn(j,k) = gauss(X(k),aInfinal(j,k),cInfinal(j,k));    % Better results
                    if muIn(j,k)<tol
                       muIn(j,k) = tol;
                    end
                    if muIn(j,k)>1.0-tol
                       muIn(j,k) = 1.0-tol;
                    end
                end            
            end
            %Layer 2: Calculation of w 
            for j=1:nRules
                    w(j) = 1.0;
                for k=1:nInputs
                    w(j) = w(j).*muIn(indexTable(j,k),k);
                end
            end
            % Layer 3: Normalization
            sumW = sum(w);
            if abs(sumW)<tol
                sumW=tol;
            end
            for j=1:nRules
               wn(j) = w(j)/sumW;
            end
            % Layer 4: Cálculo de wn_i*f_i 
            for j=1:nRules
               fi(j)   = invSigmoid(w(j),aOutfinal(j,1),cOutfinal(j,1)); % Better results
               fi_wn(j) = fi(j)*wn(j); 
            end
            % Layer 5: Final layer
             f = sum(fi_wn);
             Yfree(z+1) = f;
     	end
          
        % ANFIS (Ystep)
        for YSTEP=1:1
            % ukk1
            for UKK1=1:1
              if ukk(1,i-1)>90
                  ddu(1) = +20;
              else
                  ddu(1) = -20;
              end
              
              % ANFIS
              for z=1:Np
                % Prelayer: Input variables
                if z==1
                    X(4) = xkk(1,i-1);
                else
                    X(4) = Ystep(z-1,1);
                end
                X(1) = ukk(1,i-1)+ddu(1);
                X(2) = ukk(2,i-1);
                X(3) = ukk(3,i-1);               
                % Layer 1: Fuzzification 
                for k=1:nInputs
                    for j=1:nFuzzy
                        muIn(j,k) = gauss(X(k),aInfinal(j,k),cInfinal(j,k));    % Better results
                        if muIn(j,k)<tol
                           muIn(j,k) = tol;
                        end
                        if muIn(j,k)>1.0-tol
                           muIn(j,k) = 1.0-tol;
                        end
                    end            
                end
                %Layer 2: Calculation of w 
                for j=1:nRules
                        w(j) = 1.0;
                    for k=1:nInputs
                        w(j) = w(j).*muIn(indexTable(j,k),k);
                    end
                end
                % Layer 3: Normalization
                sumW = sum(w);
                if abs(sumW)<tol
                    sumW=tol;
                end
                for j=1:nRules
                   wn(j) = w(j)/sumW;
                end
                % Layer 4: Cálculo de wn_i*f_i 
                for j=1:nRules
                   fi(j)   = invSigmoid(w(j),aOutfinal(j,1),cOutfinal(j,1)); % Better results
                   fi_wn(j) = fi(j)*wn(j); 
                end
                % Layer 5: Final layer
                 f = sum(fi_wn);
                 Ystep(z,1) = f;
              end
            end
            
            % ukk2
            for UKK2=1:1
              if ukk(2,i-1)>50
                  ddu(2) = -25;
              else
                  ddu(2) = +25;
              end
              % ANFIS
              for z=1:Np
                % Prelayer: Input variables
                if z==1
                    X(4) = xkk(1,i-1);
                else
                    X(4) = Ystep(z-1,2);
                end
                X(1) = ukk(1,i-1);
                X(2) = ukk(2,i-1)+ddu(2);
                X(3) = ukk(3,i-1);               
                % Layer 1: Fuzzification 
                for k=1:nInputs
                    for j=1:nFuzzy
                        muIn(j,k) = gauss(X(k),aInfinal(j,k),cInfinal(j,k));    % Better results
                        if muIn(j,k)<tol
                           muIn(j,k) = tol;
                        end
                        if muIn(j,k)>1.0-tol
                           muIn(j,k) = 1.0-tol;
                        end
                    end            
                end
                %Layer 2: Calculation of w 
                for j=1:nRules
                        w(j) = 1.0;
                    for k=1:nInputs
                        w(j) = w(j).*muIn(indexTable(j,k),k);
                    end
                end
                % Layer 3: Normalization
                sumW = sum(w);
                if abs(sumW)<tol
                    sumW=tol;
                end
                for j=1:nRules
                   wn(j) = w(j)/sumW;
                end
                % Layer 4: Cálculo de wn_i*f_i 
                for j=1:nRules
                   fi(j)   = invSigmoid(w(j),aOutfinal(j,1),cOutfinal(j,1)); % Better results
                   fi_wn(j) = fi(j)*wn(j); 
                end
                % Layer 5: Final layer
                 f = sum(fi_wn);
                 Ystep(z,2) = f;
              end
            end
            
            % ukk3
            for UKK3=1:1
              if ukk(3,i-1)>50
                  ddu(3) = -35;
              else
                  ddu(3) = +35;
              end
              % ANFIS
              for z=1:Np
                % Prelayer: Input variables
                if z==1
                    X(4) = xkk(1,i-1);
                else
                    X(4) = Ystep(z-1,3);
                end
                X(1) = ukk(1,i-1);
                X(2) = ukk(2,i-1);
                X(3) = ukk(3,i-1)+ddu(3);                
                % Layer 1: Fuzzification 
                for k=1:nInputs
                    for j=1:nFuzzy
                        muIn(j,k) = gauss(X(k),aInfinal(j,k),cInfinal(j,k));    % Better results
                        if muIn(j,k)<tol
                           muIn(j,k) = tol;
                        end
                        if muIn(j,k)>1.0-tol
                           muIn(j,k) = 1.0-tol;
                        end
                    end            
                end
                %Layer 2: Calculation of w 
                for j=1:nRules
                        w(j) = 1.0;
                    for k=1:nInputs
                        w(j) = w(j).*muIn(indexTable(j,k),k);
                    end
                end
                % Layer 3: Normalization
                sumW = sum(w);
                if abs(sumW)<tol
                    sumW=tol;
                end
                for j=1:nRules
                   wn(j) = w(j)/sumW;
                end
                % Layer 4: Cálculo de wn_i*f_i 
                for j=1:nRules
                   fi(j)   = invSigmoid(w(j),aOutfinal(j,1),cOutfinal(j,1)); % Better results
                   fi_wn(j) = fi(j)*wn(j); 
                end
                % Layer 5: Final layer
                 f = sum(fi_wn);
                 Ystep(z,3) = f;
              end
            end
        end
        
        % Gij
        for nn=1:mInputs
        	for j=1:Np
                gi(j) = (Ystep(j,nn)-Yfree(j))/ddu(nn);
                %gi(j) = (Ystep(j,nn)-xkk(1,i))/ddu(nn);
            end
            for j=1:Np
                for k=1:Nc
                    if k>j
                        Gij(j,k+(nn-1)*Nc) = 0;
                    else
                        Gij(j,k+(nn-1)*Nc) = gi(j-k+1);
                    end
                end
            end 
        end
        
        % d) Optimal control sequence calculation
        duk = (((Gij'*Q*Gij+R)^(-1))*Gij'*Q*(Wref-Yfree));
        
        % e) Determination of ut
        % ukk (fuzzy)
        ukk(1,i) = (ukk(1,i-1)+006.000*duk(1,1));
        ukk(2,i) = (ukk(2,i-1)+001.000*duk(1+Nc,1));
        ukk(3,i) = (ukk(3,i-1)+002.000*duk(1+2*Nc,1));
      
        % Limit
        for LIMIT=1:1
            if ukk(1,i) < minOut
                ukk(1,i) = minOut;
            end
            if ukk(2,i) < minOut
                ukk(2,i) = minOut;
            end
            if ukk(3,i) < minOut
                ukk(3,i) = minOut;
            end

            if ukk(1,i) > maxOut
                ukk(1,i) = maxOut;
            end
            if ukk(2,i) > maxOut
                ukk(2,i) = maxOut;
            end
            if ukk(3,i) > maxOut
                ukk(3,i) = maxOut;
            end
        end
        
      end
      
      % Scaling to Force (N)
      UTT(1,i) = Mapping(ukk(1,i),minOut,maxOut,minUT1,maxUT1); % @(Input,minIn,maxIn,minOut,maxOut)
      UTT(2,i) = Mapping(ukk(2,i),minOut,maxOut,minUT2,maxUT2); % ((Input-minIn)./(maxIn-minIn)).*(maxOut-minOut)+minOut;
      UTT(3,i) = Mapping(ukk(3,i),minOut,maxOut,minUT3,maxUT3); %
  end

% Step 11: Control results
fprintf('\n\nStep 11: Control results...');
figure(10)
cla
hold on
plot(ttt, ukk(1,:),'k');
plot(ttt, ukk(2,:),'-.k');
plot(ttt, ukk(3,:),'--k');
plot(ttt, wref);
plot(ttt, xkk(1,:));
plot(ttt, yee);
xlim([3 nctrlData-Np].*step)
grid minor
title('Control results'); xlabel('Time (s)'); ylabel('% Fuzzy');
legend('Control action: ukk1','Control action: ukk2','Control action: ukk3','Reference: wref','Process output: xkk','Model output: yee');
        
%% Step 12: Full system behavior
fprintf('\n\nStep 12: Full system behavior...');
lineWidth = 1.5;

% a) Figure(1): Inputs
    figure(11)
    subplot(5,2,1);
    cla
    grid minor
    hold on;
    plot(ttt,UU(1,:),'-','LineWidth',lineWidth);
    legend('u1');
    title('Starboard motor force Tds'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,2,3);
    cla
    grid minor
    hold on;
    plot(ttt,UU(2,:),'-','LineWidth',lineWidth);
    legend('u2');
    title('Port motor force Tsp'); %xlabel('Time (s)'); 
    ylabel('(N)');
    
    subplot(5,2,2);
    cla
    grid minor
    hold on;
    plot(ttt,UTT(1,:),'-','LineWidth',lineWidth);
    legend('uT1');
    title('Motor max thrust Tmax'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,2,4);
    cla
    grid minor
    hold on;
    plot(ttt,UTT(2,:),'-','LineWidth',lineWidth);
    legend('uT2');
    title('Motor thrust difference ALPHA'); %xlabel('Time (s)'); 
    ylabel('[]');

    subplot(5,2,[5 6]);
    cla
    grid minor
    hold on;
    plot(ttt,UTT(3,:).*RadToDeg,'-','LineWidth',lineWidth);
    legend('u3');
    title('Propeller angle \mu'); %xlabel('Time (s)');  
    ylabel('(°)');

    subplot(5,2,[7 8]);
    cla
    grid minor
    hold on;
    plot(ttt,Vtot,'-','LineWidth',lineWidth);
    %legend('');
    title('Net speed'); xlabel('Time (s)');  ylabel('(m/s)');

    subplot(5,2,[9 10]);
    cla
    grid minor
    hold on;
    plot(ttt,wTot.*RadToDeg,'-','LineWidth',lineWidth);
    %legend('');
    title('Net angular speed'); xlabel('Time (s)');  ylabel('(°/s)');

% b) Figure(2): Velocities
    figure(12)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(ttt,XX(1,:),'-','LineWidth',lineWidth);
    title('u = X axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(ttt,XX(2,:),'-','LineWidth',lineWidth);
    title('v = Y axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(ttt,XX(3,:),'-','LineWidth',lineWidth);
    title('w = Z axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(ttt,XX(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('p = X axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(ttt,XX(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('q = Y axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(ttt,XX(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('r = Z axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
% c) Figure(3): Position and orientation
    figure(13)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(ttt,YY(1,:),'-','LineWidth',lineWidth);
    title('x = X axis Position'); xlabel('Time (s)'); ylabel('Position (m)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(ttt,YY(2,:),'-','LineWidth',lineWidth);
    title('y = Y axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(ttt,YY(3,:),'-','LineWidth',lineWidth);
    title('z = Z axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(ttt,YY(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\phi = X axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(ttt,YY(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\theta = Y axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(ttt,YY(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\psi = Z axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');

% d) Figure(4): Forces
    figure(14)
    subplot(2,2,1);
    cla
    grid minor
    hold on;
    plot(ttt,Fd(1,:),'-','LineWidth',lineWidth);
    plot(ttt,Fd(2,:),'-','LineWidth',lineWidth);
    plot(ttt,Fd(3,:),'-','LineWidth',lineWidth);
    plot(ttt,Fd(4,:),'-','LineWidth',lineWidth);
    plot(ttt,Fd(5,:),'-','LineWidth',lineWidth);
    plot(ttt,Fd(6,:),'-','LineWidth',lineWidth);
    title('Fd(X(t))'); xlabel('Time (s)'); ylabel('Dynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,2);
    cla
    grid minor
    hold on;
    plot(ttt,G(1,:),'-','LineWidth',lineWidth);
    plot(ttt,G(2,:),'-','LineWidth',lineWidth);
    plot(ttt,G(3,:),'-','LineWidth',lineWidth);
    plot(ttt,G(4,:),'-','LineWidth',lineWidth);
    plot(ttt,G(5,:),'-','LineWidth',lineWidth);
    plot(ttt,G(6,:),'-','LineWidth',lineWidth);
    title('G(X(t))'); xlabel('Time (s)'); ylabel('Gravity force (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,3);
    cla
    grid minor
    hold on;
    plot(ttt,A(1,:),'-','LineWidth',lineWidth);
    plot(ttt,A(2,:),'-','LineWidth',lineWidth);
    plot(ttt,A(3,:),'-','LineWidth',lineWidth);
    plot(ttt,A(4,:),'-','LineWidth',lineWidth);
    plot(ttt,A(5,:),'-','LineWidth',lineWidth);
    plot(ttt,A(6,:),'-','LineWidth',lineWidth);
    title('A(X(t))'); xlabel('Time (s)'); ylabel('Aerodynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,4);
    cla
    grid minor
    hold on;
    plot(ttt,P(1,:),'-','LineWidth',lineWidth);
    plot(ttt,P(2,:),'-','LineWidth',lineWidth);
    plot(ttt,P(3,:),'-','LineWidth',lineWidth);
    plot(ttt,P(4,:),'-','LineWidth',lineWidth);
    plot(ttt,P(5,:),'-','LineWidth',lineWidth);
    plot(ttt,P(6,:),'-','LineWidth',lineWidth);
    title('P(u(t))'); xlabel('Time (s)'); ylabel('Propulsion forces (N)');
    legend('x','y','z','roll','pitch','yaw');  
%%              
fprintf('\nEND\n');
toc