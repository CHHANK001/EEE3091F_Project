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
Vp = double(Vline / sqrt(3));


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
Ns = (120*f)/p;                     %Ns is defined twice need to fix that
    %disp(['Ns  = ' num2str(Ns)])
ws = Ns*(2*pi/60);

% Question 2 
% Cerate s matrix
Smin = 1;
Sstep = 1;
Smax = 100;
s = Smin:Sstep:Smax;
% create Ns matrix
Ns = flip(s);           % this is wrong and needs to be changed


    
% Calc T
Tm1 = 3*(1/ws).*((Vth(1)^2) ./ ((Rth(1)+R2(1) ./ (0.01*Ns).^2)+((Xth(1) + X2(1)).^2))).*(R2(1)./(0.01*Ns));
Tm2 = 3*(1/ws).*((Vth(2)^2) ./ ((Rth(2)+R2(2) ./ (0.01*Ns).^2)+((Xth(2) + X2(2)).^2))).*(R2(2)./(0.01*Ns));

% create a figure window with subplots
figure;
subplot(3,1,1);
% Plot Torque vs speed
plot(s,Tm1, 'b-',s,Tm2,'r-')                    % plots are showing %Ns but need to be showing n 
title('Torque vs Speed Characteristic')
xlabel('%Ns')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor')

% Question 3 Stator current vs Speed
Z11 = R1(1) + complex(0, X1(1)) + (complex(0,Xm(1)).*((R2(1)./(0.01*Ns))+complex(0,X2(1)))) ./ ((R2(1)./(0.01*Ns)) + complex(0, Xm(1)+X2(1))); %eqn is correct but outputting wrong values
Z12 = R1(2) + complex(0, X1(2)) + (complex(0,Xm(2)).*((R2(2)./(0.01*Ns))+complex(0,X2(2)))) ./ ((R2(2)./(0.01*Ns)) + complex(0, Xm(2)+X2(2)));

disp(['Z11 = ' num2str(Z11)])
%disp(['Z12 = ' num2str(Z12)])
I11 = Vp ./ abs(Z11);
I12 = Vp ./ abs(Z12);

% Plot current vs speed characteristic
subplot(3,1,2)
plot(s,I11, 'b-',s,I12,'r-')
title('Stator current vs Speed Characteristic')
xlabel('%s')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor')

% Question 4
theta1 = atan(imag(Z11)./real(Z11));
pf1 = cos(theta1);
theta2 = atan(imag(Z12) ./ real(Z12));
pf2 = cos(theta2);

%plot pf vs %Ns
subplot(3,1,3)
plot(s,pf1, 'b-',s,pf2,'r-')
title('Power Factor vs Speed Characteristic')
xlabel('%Ns')
ylabel('Power Factor')
legend('SE Motor', 'EE Motor')

% Question 5





