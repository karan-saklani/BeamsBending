function [x,y] = center(n)
A = zeros(n);
xi = zeros(n);
yi = zeros(n);
for i = 1:n
    A(i) = input('Enter the Area ');
    xi(i) = input('Enter x distance from center ');
    yi(i) = input('Enter y distance from center ');
end
x = sum(A.*xi)/sum(A);
y = sum(A.*yi)/sum(A);
end

function Ix = MIx(n)
Ix = 0;
for i = 1:n
opt = input('Enter Shape ','s');
if opt == "r"
    [I,A] = I_Rect('x');
    y = input('Enter y distance of centroid from center ');
    Ix = Ix + I + A*y^2;
elseif opt == "sc"
    [I,A] = I_SC('x');
    y = input('Enter y distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ' );
    Ix = Ix + I + A*y^2;
elseif opt == "c"
    [I,A] = I_Circle('x');
    y = input('Enter y distance of centroid from center ');
    Ix = Ix + I + A*y^2;
elseif opt == "rt"
    [I,A] = I_RT('x');
    y = input('Enter y distance of centroid from center ( which is at height of h/3 from base ');
    Ix = Ix + I + A*y^2;
end
end
end

function Iy = MIy(n)
Iy = 0;
for i = 1:n
opt = input('Enter Shape','s');
if opt == "r"
    [I,A] = I_Rect('y');
    x = input('Enter x distance of centroid from center ');
    Iy = Iy + I + A*x^2;
elseif opt == "sc"
    [I,A] = I_SC('y');
    x = input('Enter x distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ' );
    Iy = Iy + I + A*x^2;
elseif opt == "c"
    [I,A] = I_Circle('y');
    x = input('Enter x distance of centroid from center ');
    Iy = Iy + I + A*x^2;
elseif opt == "rt"
    [I,A] = I_RT('y');
    x = input('Enter x distance of centroid from center ( which is at distance of b/3 from right angled vortex) ');
   
    Iy = Iy + I + A*x^2;
end
end
end

function Ixy = MIxy(n)
Ixy = 0;
for i = 1:n
opt = input('Enter Shape','s');
if opt == "r"
    [I,A] = I_Rect('xy');
    x = input('Enter x distance of centroid from center ');
    y = input('Enter y distance of centroid from center ');
    Ixy = Ixy + I + A*x*y;
elseif opt == "sc"
    [I,A] = I_SC('xy');
    x = input('Enter x distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ' );
    y = input('Enter y distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ' );
    Ixy = Ixy + I + A*x*y;
elseif opt == "c"
    [I,A] = I_Circle('xy');
    x = input('Enter x distance of centroid from center ');
    y = input('Enter y distance of centroid from center ')
    Ixy = Ixy + I + A*x*y;
elseif opt == "rt"
    [I,A] = I_RT('xy');
    x = input('Enter x distance of centroid from center ( which is at distance of b/3 from right angled vortex) ');
    y = input('Enter y distance of centroid from center ( which is at height of h/3 from base ');
    Ixy = Ixy + I + A*x*y;
end
end
end

function [Ix,Iy,Ixy] = MI(n)
Ix = 0;
Iy = 0;
Ixy = 0;
for i = 1:n
    opt = input('Enter Shape [Rectangle(r),SemiCircle(sc)] ','s');
    if opt == "r"
        h = input('Enter value of h distance along y axis) ');
        b = input('Enter value of b distance along x axis');
        [I1,A] = I_Rect_given('x',h,b);
        y = input('Enter y distance of centroid from center');
        Ix = Ix + I1 + A*y^2;
        [I2,A] = I_Rect_given('y',h,b);
        x = input('Enter x distance of centroid from center');
        Iy = Iy + I2 + A*x^2;
        [I3,A] = I_Rect_given('xy',h,b);
        Ixy = Ixy + I3 + A*x*y;
    elseif opt == "sc"
        R = input('Enter Value of Radius');
        [I1,A] = I_SC_given('x',R);
        y = input('Enter y distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ');
        Ix = Ix + I1 + A*y^2;
        [I2,A] = I_SC_given('y');
        x = input('Enter x distance of centroid from center (which is 4R/3pi along the height of the semi-circle center) from center ');
        Iy = Iy + I2 + A*x^2;
        [I3,A] = I_SC_given('xy');
        Ixy = Ixy + I3 + A*x*y;
    elseif opt == "rt"
        h = input('Enter value of h distance along y axis) ');
        b = input('Enter value of b distance along x axis');
        [I1,A] = I_RT_given('x',b,h);
        y = input('Enter y distance of centroid from center ( which is at height of h/3 from base ');
        Ix = Ix + I1 + A*y^2;
        [I2,A] = I_RT_given('y',b,h);
        x = input('Enter x distance of centroid from center ( which is at distance of b/3 from hypotenuese vortex) ');
        Iy = Iy + I2 + A*x^2;
        [I3,A] = I_RT_given('xy',b,h);
        Ixy = Ixy + I3 + A*x*y;
    elseif opt == "c"
        R = input('Enter Value of Radius');
        [I1,A] = I_Circle_given('x',R);
        y = input('Enter y distance of centroid from center ');
        Ix = Ix + I1 + A*y^2;
        [I2,A] = I_Circle_given('y',R);
        x = input('Enter x distance of centroid from center ');
        Iy = Iy + I2 + A*x^2;
        [I3,A] = I_Circle_R('xy',R);
        Ixy = Ixy + I3 + A*x*y;        
    end
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

function [I,A] = I_SC(dir)
R = input('Enter value of Radius');
if dir == 'x'
    I = pi*R^4*(1/8 - 8/(9*pi^2));
elseif dir == 'y'
    I = pi*R^4/8;
elseif dir == 'r'
    I = pi*R^4*(1/4 - 8/(9*pi^2));
elseif dir == "xy"
    I = 0;
end
A = (pi*R^2)/2;
end

function [I,A] = I_RT(dir)
h = input('Enter value of h distance along y axis) ');
b = input('Enter value of b distance along x axis');
if dir == 'x'
    I = b*h^3/36;
elseif dir == 'y'
    I = h*b^3/36;
elseif dir == 'r'
    I = (b*h^3 + h*b^3)/36;
elseif dir == "xy"
    I = -(b^2*h^2)/72;
end
A = 0.5*b*h;
end

function [I,A] = I_Circle(dir)
R = input('Enter value of Radius');
if dir == 'x'
    I = pi*R^4/4;
elseif dir == 'y'
    I = pi*R^4/4;
elseif dir == 'r'
    I = pi*R^4/2;
elseif dir == "xy"
    I = 0;
end
A = (pi*R^2);
end

% Given functions
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

function [I,A] = I_SC_given(dir,R)
R = input('Enter value of Radius');
if dir == 'x'
    I = pi*R^4*(1/8 - 8/(9*pi^2));
elseif dir == 'y'
    I = pi*R^4/8;
elseif dir == 'r'
    I = pi*R^4*(1/4 - 8/(9*pi^2));
elseif dir == "xy"
    I = 0;
end
A = (pi*r^2)/2;
end

function [I,A] = I_RT_given(dir,b,h)
if dir == 'x'
    I = b*h^3/36;
elseif dir == 'y'
    I = h*b^3/36;
elseif dir == 'r'
    I = (b*h^3 + h*b^3)/36;
elseif dir == "xy"
    I = -(b^2*h^2)/72;
end
A = 0.5*b*h;
end

function [I,A] = I_Circle_given(dir,R)
if dir == 'x'
    I = pi*R^4/4;
elseif dir == 'y'
    I = pi*R^4/4;
elseif dir == 'r'
    I = pi*R^4/2;
elseif dir == "xy"
    I = 0;
end
A = (pi*R^2);
end
