clear;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Energy Efficient (EE)Motor parameters')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
 
f=50;               %Supply frequency [Hz]
p=4;                %Number of poles
V1=380/sqrt(3);     %Supply voltage [phase]
R1=1.5;             %Stator winding resistance [ohms/phase]
X1=3.642;           %Stator winding leakage reactance [ohms/phase]
Xm=72.252;          %Stator winding magnetising reactance [ohms/phase]
X2p=3.642;          %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2p=1.994;          %Rotor winding resistance reffered to stator [ohms/phase]

fprintf('f=%f\n',f);
fprintf('p=%f\n',p);
fprintf('V1=%f\n',V1);
fprintf('R1=%f\n',R1);
fprintf('X1=%f\n',X1);
fprintf('Xm=%f\n',Xm);
fprintf('X2p=%f\n',X2p);
fprintf('R2p=%f\n',R2p);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 1:Thevenin Equiv Cct Parameters for EE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Vth=Xm/sqrt(R1^2+(X1+Xm)^2)*V1;         %Thevenin equiv voltage source [V] (Equ 5.45 - Sen)
Zth=1i*Xm*(R1+1i*X1)/(R1+1i*(X1+Xm));   %Thevenin equiv impedance
Rth=real(Zth);                          %Thevenin equiv resistance [ohms]
Xth=imag(Zth);                          %Thevenin equiv reactance [ohms]

fprintf('Vth=%f\n',Vth);
fprintf('Rth=%f\n',Rth);
fprintf('Xth=%f\n',Xth);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('QUESTION 2:Torque versus speed characteristics for EE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

ns=120*f/p;         %Synchronous speed [rpm]

ws=2*pi*ns/60;      %Synchronous speed [rad/sec]
s=0.0005:0.0005:1;  %Slip [pu]
disp(["s = " num2str(s)]);

n=(1-s)*ns;         %Rotor speed [rpm]
disp(["n = " num2str(n)]);

w=2*pi*n/60;        %Rotor speed [rad/sec]
 
Tmech=3/ws*Vth^2./((Rth+R2p./s).^2+(Xth+X2p)^2).*R2p./s;    %Total Tmech = {3*(Equ5.54 - Sen)}
disp(["Tmech = " num2str(Tmech)]);

%Plot the Torque vs speed characteristics of motor:
subplot(2,2,1),
plot(n,Tmech),xlabel('n [rpm]'),ylabel('Torque [Nm]'),...
   title('Torque vs speed'),grid on,...
   hold on