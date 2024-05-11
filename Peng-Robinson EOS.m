% PENG ROBINSON EQUATION OF STATE
% DATE : 10-01-21
% NAME : SANTOSH M. ZOL


clc
format short

T = input('enter the value of temperature ='); % temp in degree celcius
T_k = T + 273.15; % temp in degree kelvin
p = 100000; % pressure in pa

% component1 = methane; component2 = ethane; component3 = propane
Y = [0.5 0.3 0.2];
M = [16.04 30.7 44.1];
Tc = [190.6 305.4 370.0];
Pc = [4600000 4874000 4244000]; % pressure in pa
R = 8.314; 
omega = [0.008 0.099 0.152];
n = length(Tc);

a = zeros(1,n);
b = zeros(1,n);
k = zeros(1,n);
Tr = zeros(1,n);
alpha = zeros(1,n);

for i=1:n
    a(1,i) = 0.45724*((R^2)*(Tc(1,i))^2)/Pc(1,i);
    b(1,i) = 0.07780*(R*Tc(1,i))/Pc(1,i);
    k(1,i) = 0.37464 + 1.54226*(omega(1,i)) - 0.26992*(omega(1,i))^2;
    Tr(1,i) = T_k/(Tc(1,i));
    alpha(1,i) = (1+k(1,i)*(1-sqrt(Tr(1,i))))^2;
    
end

a_alpha = zeros(n);
for i=1:n
    for j=1:n
        a_alpha(i,j) = sqrt(a(1,i)*a(1,j)*alpha(1,i)*alpha(1,j));
    end
end


% claculating a_mix
a_mix = 0;
for i=1:n
    for j=1:n
        a_mix = a_mix + Y(1,i)*Y(1,j)*a_alpha(i,j);
    end
end

%calculating b_mix
b_mix = 0;
for i=1:n
    b_mix = b_mix + Y(1,i)*b(1,i);
end

% BY REARRANGENMENT WE GET THE PENG ROBINSON EQUATION IN POLYNOMIAL FORM
% V^3 + A*V^2 + B*V + D = 0
% calculating values of A, B and D

 A = b_mix - (R*T_k)/p;
 B = ((a_mix - 2*R*T_k*b_mix)/p) - 3*(b_mix)^2;
 D = (b_mix)^3 + ((R*T_k)*(b_mix)^2)/p - (a_mix*b_mix)/p;
 
 % finding the root of polynomial
 
h = 0.01; % Initial guess
t = 0.001; % Enter the tolerance value
y = @(v) v.^3 + A*v.^2 + B*v + D;
z = @(v) 3*v.^2 + 2*A*v + B;
if y(h) == 0
    fprintf('The molar volume of gaseous mixture is %f\n', h);
end
while abs(y(h))>0
    c = h - (y(h)/z(h));
    Mol_D = 0.001*(sum(M.*Y))/c;
    
    if abs(h-c)<t
        fprintf('Molar Volume of gaseous mixture is %f\n',c)
        fprintf('molar density of gaseous mixture is %f\n',Mol_D)
        break
    end
    h=c;
    
end