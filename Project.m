clc                                             %Clears the command window
% EEE3091 Project 
% Authors:  Ankush Chohan (CHHANK001)
%           Ashik John (JHNASH009)
% Date of last revision: 31/03/2023


% Section A
% Define variables for both motors in an array [SE,EE]

R1 = double([2.087, 1.500]);                    %Stator winding resistance [ohms/phase]
X1 = double([4.274, 3.642]);                    %Stator winding leakage reactance [ohms/phase]
Xm = double([66.560, 72.252]);                  %Stator winding magnetising reactance [ohms/phase]
X2 = double([4.274, 3.642]);                    %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2 = double([2.122, 1.994]);                    %Rotor winding resistance reffered to stator [ohms/phase]
Prot = double([134.669, 88.924]);
Vline = 380;
f = 50;                                         %Supply frequency [Hz]
p = 4;                                          %Number of poles
Vp = double(Vline / sqrt(3));                   %Supply voltage [phase]


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 1:Thevenin Equiv Cct Parameters for both motors [SE,EE]:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


Vth = (Xm ./ sqrt(R1.^2 + (X1+Xm).^2)).*Vp;     %Thevenin equiv voltage source [V]


Voc = (complex(0,Xm)./(R1 + complex(0,X1+Xm))).*Vp;
%disp(['Voc = ' num2str(Voc)])

Isc = Vp./(R1 + complex(0, X1));
%disp(['Isc = ' num2str(Isc)])
  
Zth = Voc ./ Isc;                               %Thevenin equiv impedance
%disp(['Zth = ' num2str(Zth)])

Rth = real(Zth);                                %Thevenin equiv resistance [ohms]
Xth = imag(Zth);                                %Thevenin equiv reactance [ohms]

disp(['Vth = ' num2str(Vth)]);
disp(['Rth = ' num2str(Rth)]);
disp(['Xth = ' num2str(Xth)]);


% Question 2
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 2:Torque vs Speed characteristics for both motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


Ns = (120*f)/p;                                 %Synchronous speed [rpm]
ws = Ns*(2*pi/60);                              %Synchronous speed [rad/sec]


% Create s matrix
s = 0.0005:0.0005:1;                            %Slip [pu]
% create Ns matrix   
n = (1-s).*Ns;                                  %Rotor speed [rpm]
w = n.*(2*pi/60);                               %Rotor speed [rad/sec]

    
% Calc T
Tm1 = 3*(1/ws).*((Vth(1)^2) ./ (((Rth(1)+(R2(1) ./ (s))).^2)+((Xth(1) + X2(1)).^2))).*(R2(1)./(s));
Tm2 = 3*(1/ws).*((Vth(2)^2) ./ (((Rth(2)+(R2(2) ./ (s))).^2)+((Xth(2) + X2(2)).^2))).*(R2(2)./(s));

figure(1);
% Plot Torque vs speed
plot(n,Tm1, 'b-',n,Tm2,'r-')                   
title('Torque vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor')



disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 3:Stator current vs Speed characteristics for both motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Z11 = R1(1) + complex(0, X1(1)) + (complex(0,Xm(1)).*((R2(1)./(s))+complex(0,X2(1)))) ./ ((R2(1)./(s)) + complex(0, Xm(1)+X2(1))); 
Z12 = R1(2) + complex(0, X1(2)) + (complex(0,Xm(2)).*((R2(2)./(s))+complex(0,X2(2)))) ./ ((R2(2)./(s)) + complex(0, Xm(2)+X2(2)));

I11 = Vp ./ abs(Z11);                           %Current in motor 1
I12 = Vp ./ abs(Z12);                           %Current in motor 2

% Plot current vs speed characteristic
figure(2);
plot(n,I11, 'b-',n,I12,'r-')
title('Stator current vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Stator current(A)')
legend('SE Motor', 'EE Motor')


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 4:Power Factor vs Speed characteristics for both motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

theta1 = atan(imag(Z11)./real(Z11));
pf1 = cos(theta1);
theta2 = atan(imag(Z12) ./ real(Z12));
pf2 = cos(theta2);

%plot pf vs speed
figure(3);
%subplot(3,1,3)
plot(n,pf1, 'b-',n,pf2,'r-')
title('Power Factor vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Power Factor')
legend('SE Motor', 'EE Motor')


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 5:Power vs Speed characteristics for both motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% Question 5 rotational loss excluded
% Motor 1                                       % Eqns from slide 21 Week 6
Pin1 = 3*Vp*I11.*pf1;                           % Input power
P1cu1 = 3*(I11.^2)*R1(1);                       % Stator copper loss
P2cu1 = 3*(I11.^2)*R2(1);                       % Rotor copper loss
Pag1 = Pin1 - P1cu1;                            % Airgap power
Pshaft1 = Pag1 - P2cu1;                         % Shaft power


% Motor 2
Pin2 = 3*Vp.*I12 .* pf2;                        % Input power
P1cu2 = 3.*(I12.^2)*R1(2);                      % Stator copper loss
P2cu2 = 3.*(I12.^2)*R2(2);                      % Rotor copper loss
Pag2 = Pin2 - P1cu2;                            % Airgap power
Pshaft2 = Pag2 - P2cu2;                         % Shaft power

% Create figure 4 
figure(4);
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

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 6: Efficiency vs Speed characteristics for both motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

% Motor 1
n1 = Pshaft1./Pin1;                             % Calculate efficiency

% Motor 2
n2 = Pshaft2./Pin2;                             % Calculate efficiency

% Plot efficiency vs speed for motor 1
% Create figure 5
figure(5);
% Plot efficiency for motor 1
plot(n,n1,n,n2);
title('Efficiency vs Speed Characteristic');
xlabel('Speed [rpm]');
ylabel('Eff');
legend('SE Motor','EE Motor');





disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 7: Torque vs Speed characteristics for both motors + Centrifugal Pump load:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% Question 7

k = 946.88 * 10^-6;                             % Define constant variable
Tlo = k.*(w.^2);                                % Calculate load torque

% Create figure 6
figure(6)

% Plot Torque vs speed with load torque displayed aswell 
plot(n,Tm1,n,Tm2,n, Tlo)                   
title('Torque vs Speed Characteristic')
xlabel('Speed [rpm]')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor','Tlo')



