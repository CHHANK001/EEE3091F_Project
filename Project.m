clc %Clears the command window
% EEE3091 Project 
% Authors:  Ankush Chohan (CHHANK001)
%           Ashik John (JHNASH009)
% Date of last revision: 30/03/2023

% Section A Q1

% Define variables
R1 = double([2.087, 1.500]);
X1 = double([4.274, 3.462]);
Xm = double([66.560, 72.252]);
X2 = double([4.274, 3.642]);
R2 = double([2.122, 1.994]);
Prot = double([134.669, 88.924]);
Vline = 380;
f = 50;
p = 3;
Vp = double(Vline / sqrt(3));       % Calculates phase voltage


% Calculate Vth
Vth = (Xm ./ sqrt(R1.^2 + (X1+Xm).^2)).*Vp;
disp(['Vth ' num2str(Vth)])


% Calculate Voc
Voc = complex(0,Xm)./(R1 + complex(0,X1+Xm)).*Vp;
disp(['Voc = ' num2str(Voc)])

% Calculate Isc
Isc = Vp./(R1 + complex(0, X1));
disp(['Isc = ' num2str(Isc)])
  
% Calculate Zth
Zth = Voc ./ Isc;
disp(['Zth = ' num2str(Zth)])

% Separate out Rth and Zth 
Rth = real(Zth);
Xth = imag(Zth);
disp(['Rth = ' num2str(Rth)])
disp(['Xth = ' num2str(Xth)])

% Question 2

% find Ns and ws
Ns = (120*f)/p;                    
ws = Ns*(2*pi/60);          % Ns in rad/s

% Question 2 
% Cerate s matrix
Smin = 1;
Sstep = 1;
Smax = 100;
a = Smin:Sstep:Smax;
s = 0.01*a;
% create Ns matrix
f = flip(s);          
n = s.*Ns;
w = n.*(2*pi/60);
disp(['w = ' num2str(w)])
    
% Calc T
Tm1 = 3*(1/ws).*((Vth(1)^2) ./ ((Rth(1)+R2(1) ./ (f).^2)+((Xth(1) + X2(1)).^2))).*(R2(1)./(f));
Tm2 = 3*(1/ws).*((Vth(2)^2) ./ ((Rth(2)+R2(2) ./ (f).^2)+((Xth(2) + X2(2)).^2))).*(R2(2)./(f));

% Create figure 1 window with subplots
figure(1);
subplot(3,1,1);
% Plot Torque vs speed
plot(n,Tm1, 'b-',n,Tm2,'r-')                   
title('Torque vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor')

% Question 3 Stator current vs Speed
Z11 = R1(1) + complex(0, X1(1)) + (complex(0,Xm(1)).*((R2(1)./(f))+complex(0,X2(1)))) ./ ((R2(1)./(f)) + complex(0, Xm(1)+X2(1))); %eqn is correct but outputting wrong values
Z12 = R1(2) + complex(0, X1(2)) + (complex(0,Xm(2)).*((R2(2)./(f))+complex(0,X2(2)))) ./ ((R2(2)./(f)) + complex(0, Xm(2)+X2(2)));

    %disp(['Z11 = ' num2str(Z11)])
    %disp(['Z12 = ' num2str(Z12)])
I11 = Vp ./ abs(Z11);           %Current in motor 1
I12 = Vp ./ abs(Z12);           %Current in motor 2

% Plot current vs speed characteristic
subplot(3,1,2)
plot(n,I11, 'b-',n,I12,'r-')
title('Stator current vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Stator current(A)')
legend('SE Motor', 'EE Motor')

% Question 4
theta1 = atan(imag(Z11)./real(Z11));
pf1 = cos(theta1);
theta2 = atan(imag(Z12) ./ real(Z12));
pf2 = cos(theta2);

%plot pf vs speed
subplot(3,1,3)
plot(n,pf1, 'b-',n,pf2,'r-')
title('Power Factor vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Power Factor')
legend('SE Motor', 'EE Motor')

% Question 5 rotational loss excluded
% Motor 1                           % Eqns from slide 21 Week 6
Pin1 = 3*Vp*I11;                    % Input power
P1cu1 = 3*(I11.^2)*R1(1);           % Stator copper loss
P2cu1 = 3*(I11.^2)*R2(1);           % Rotor copper loss
Pag1 = Pin1 - P1cu1;                % Airgap power
Pshaft1 = Pag1 - P2cu1;             % Shaft power
disp(['Pin1 = ' num2str(Pin1)])
disp(['Pshaft1 = ' num2str(Pshaft1)])

% Motor 2
Pin2 = 3*Vp.*I12;                      % Input power
P1cu2 = 3.*(I12.^2)*R1(2);             % Stator copper loss
P2cu2 = 3.*(I12.^2)*R2(2);             % Rotor copper loss
Pag2 = Pin2 - P1cu2;                   % Airgap power
Pshaft2 = Pag2 - P2cu2;                % Shaft power

% Create figure 2 
figure(2);
% Plot for Motor 1
subplot(2,1,1)
plot(n,Pin1,n,P1cu1,n, P2cu1,n, Pag1,n, Pshaft1)
title('Motor 1 Powers vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('P')
legend('Pin', 'P1cu', 'P2cu', 'Pag', 'Pshaft')

% Plot for Motor 2
subplot(2,1,2)
plot(n,Pin2,n,P1cu2,n, P2cu2,n, Pag2,n, Pshaft2)
title('Motor 2 Powers vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('P')
legend('Pin', 'P1cu', 'P2cu', 'Pag', 'Pshaft')

% Question 6

% Motor 1
n1 = Pshaft1./Pin1;        % Calculate efficiency

% Plot efficiency vs speed for motor 1
% Create figure 3
figure(3)
% Plot efficiency for motor 1
subplot(2,1,1)
plot(n,n1)
title('Motor 1 Efficiency vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Eff')

% Motor 2
n2 = Pshaft2./Pin2;         % Calculate efficiency

% Plot efficiency for motor 2
subplot(2,1,2)
plot(n,n2)
title('Motor 2 Efficiency vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Eff')

% Question 7

k = 946.88 * 10^-6;     % Define constant variable
Tlo = k.*(w.^2);        % Calculate load torque
    %disp(['Tlo = ' num2str(Tlo)])

% Create figure 4
figure(4)

% Plot Torque vs speed with load torque displayed aswell 
plot(n,Tm1,n,Tm2,n, Tlo)                   
title('Torque vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor','Tlo')



