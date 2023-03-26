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
Vth = zeros(1,2);

% Calculate Vth
for i=1.0:1.0:2.0
  
    Vth(i) = (Xm(i) / sqrt(R1(i)^2 + (X1(i)+Xm(i)^2)))*Vp;
    disp(['Vth' num2str(i) ' = ' num2str(Vth(i))])
end

% Calculate Voc
Voc = complex(zeros(1,2), zeros(1,2));
for i=1.0:1.0:2.0
    Voc(i) = ((complex(0,Xm(i)))/(R1(i) + complex(0,X1(i)+Xm(i))))*Vp;
    disp(['Voc' num2str(i) ' = ' num2str(Voc(i))])
end

% Calculate Isc
Isc = complex(zeros(1,2), zeros(1,2));
for i=1.0:1.0:2.0
    Isc(i) = Vp/(R1(i) + complex(0, X1(i)));
    disp(['Isc' num2str(i) ' = ' num2str(Isc(i))])
end

% Calculate Zth
Zth = complex(zeros(1,2), zeros(1,2));
for i=1.0:1.0:2.0
    Zth(i) = Voc(i)/Isc(i);
    disp(['Zth' num2str(i) ' = ' num2str(Zth(i))])
end

Rth = zeros(1,2);
Xth = zeros(1,2);
for i=1.0:1.0:2.0
    Rth(i) = real(Zth(i));
    Xth(i) = imag(Zth(i));
    disp(['Rth' num2str(i) ' = ' num2str(Rth(i))])
    disp(['Xth' num2str(i) ' = ' num2str(Xth(i))])
end
% Question 2

% find Ns and ws
Ns = (120*f)/p;
disp(['Ns' num2str(i) ' = ' num2str(Ns)])
ws = Ns*(2*pi/60);

% Question 2 
% Cerate s matrix
Smin = 0;
Sstep = 1;
Smax = 100;
s = Smin:Sstep:Smax;


% Extract required constants
Vth1 = Vth(1);          Vth2 = Vth(2);               
Rth1 = Rth(1);          Rth2 = Rth(2);
R21 = R2(1);            R22 = R2(2);
Xth1 = Xth(1);          Xth2 = Xth(2);
X21 = X2(1);            X22 = X2(2);
Tm1 = zeros(1,101);     
Tm2 = zeros(1,101);     
% Calc T
for i = 1.0:1.0:101
    Tm1(i) = 3*(1/ws)*((Vth1^2)/(((Rth1+R21/(0.01*s(i))^2))+((Xth1 + X21)^2)))*(R21/(0.01*s(i)));
       %a = 3*(1/ws);
    %b = (Rth1+R21/s(i))^2 + (Xth1+X21)^2;
    %c = R21 / s(i);
    %Tm1(i) = a*((Vth1^2)/b)*c;
    Tm2(i)  = 3*(1/ws)*((Vth2^2)/(((Rth2+R22/(0.01*s(i))^2))+((Xth2 + X22)^2)))*(R22/(0.01*s(i)));
end

plot(s,Tm1, 'b-',s,Tm2,'r-')
title('Torque vs Speed Characteristic')
xlabel('%slip')
ylabel('Torque(Nm)')
legend('SE Motor', 'EE Motor')




