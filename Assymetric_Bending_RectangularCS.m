phi = 0; % value of phi varies with the problem
[data,Ix,Iy,Ixy,x0,y0,alpha,theta] = MIsolve(phi)
function [data,Ix,Iy,Ixy,x0,y0,alpha,theta] = MIsolve(phi)
data = readtable("MI.xlsx");
delete MISolution.xlsx;
xi = data.xi;
yi = data.yi;
A = data.b .* data.h;
[x0,y0] = center_excel(A,xi,yi);
data.x = x0 - xi;
data.y = y0 - yi;
shape = size(data);
n = shape(1);
[Ix,Iy,Ixy] = MI_excel_all(n,data.b,data.h,data.x,data.y);
theta = atan(2*Ixy/(Iy - Ix))/2;
alpha = neutral_axis_angle(Ix,Iy,Ixy,phi);
Ix_arr = MIx_excel(n,data.y,data.h,data.b);
Iy_arr = MIy_excel(n,data.x,data.h,data.b);
Ixy_arr = MIxy_excel(n,data.x,data.y,data.h,data.b);
data.Ix = Ix_arr;
data.Iy = Iy_arr;
data.Ixy = Ixy_arr;
writetable(data,'MISolution.xlsx');
end

function [x,y] = center_excel(A,xi,yi)
x = sum(A.*xi)/sum(A);
y = sum(A.*yi)/sum(A);
end

function Ix = MIx_excel(n,y,h,b)
Ix = zeros(n,1);
for i = 1:n
[I,A] = I_Rect_given("x",h(i),b(i));
Ix(i) = I + A*y(i)^2;
end
end

function Iy = MIy_excel(n,x,h,b)
Iy = zeros(n,1);
for i = 1:n
[I,A] = I_Rect_given("y",h(i),b(i));
Iy(i) = I + A*x(i)^2;
end
end

function Ixy = MIxy_excel(n,x,y,h,b)
Ixy = zeros(n,1);
for i = 1:n
    [I,A] = I_Rect_given("xy",h(i),b(i));
    Ixy(i) =I + A*x(i)*y(i);

end
end

function [Ix,Iy,Ixy] = MI_excel_all(n,b,h,x,y)
Ix = 0;
Iy = 0;
Ixy = 0;
for i = 1:n
    [I1,A] = I_Rect_given('x',h(i),b(i));
    Ix = Ix + I1 + A*y(i)^2;
    [I2,A] = I_Rect_given('y',h(i),b(i));
    Iy = Iy + I2 + A*x(i)^2;
    [I3,A] = I_Rect_given('xy',h(i),b(i));
    Ixy = Ixy + I3 + A*x(i)*y(i);
    
end
end

function alpha = neutral_axis_angle(Ix,Iy,Ixy,phi)
tan_alpha = (Ixy - Ix*cot(phi))/(Iy - Ixy*cot(phi));
alpha = atan(tan_alpha);
end

function [theta,Ix_p,Iy_p] = principal_axis(Ix,Iy,Ixy)
disp('Note : Principal Axis : Ixy = 0');
theta = (atan(2*Ixy/(Iy - Ix)))/2;
Ix_p = (Ix + Iy)/2 + sqrt(((Ix - Iy)/2)^2 + Ixy^2);
Iy_p = (Ix + Iy)/2 - sqrt(((Ix - Iy)/2)^2 + Ixy^2);
end

function sigmaP = stress(Mx,Ix,Ixy,alpha,y,x)
sigmaP = (Mx*(y - x*tan(alpha)))/(Ix - Ixy*tan(alpha));
end

function [u,v] = deflection(P,L,phi,E,alpha,Ix,Ixy)
v = (P*L^3*sin(phi))/(3*E*(Ix - Ixy*tan(alpha)));
u = -v*tan(alpha);
end


% Useful Functions
function [I,A] = I_Rect(dir)
h = input('Enter value of h distance along y axis) ');
b = input('Enter value of b distance along x axis');
if dir == 'x'
    I = b*h^3/12;
elseif dir == 'y'
    I = h*b^3/12;
elseif dir == 'r'
    I = (b*h^3 + h*b^3)/12;
elseif dir == "xy"
    I = 0;
end
A = b*h;
end

function [I,A] = I_Rect_given(dir,h,b)
if dir == 'x'
    I = b*h^3/12;
elseif dir == 'y'
    I = h*b^3/12;
elseif dir == 'r'
    I = (b*h^3 + h*b^3)/12;
elseif dir == "xy"
    I = 0;
end
A = b*h;
end
