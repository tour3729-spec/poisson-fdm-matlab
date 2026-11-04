% Finite Difference Method for 1 x 1 m^2 for Laplace equation ( = 4 )

clear,clc

NX = 50; % Discretization on x-axis
NY = 50; % Discretization on y-axis

DX = 1/NX;
DY = 1/NY;

IMAX = NX + 1; % Maximum number of NODES on x-axis
JMAX = NY + 1; % Maximum number of NODES on y-axis

% GRID GENERATION
meshview = input('Do you want to plot grid? (type 1 if TRUE): ');

if meshview == 1
    X = 0:DX:1;
    Y = 0:DY:1;
    Y = Y';

    for i1 = 2:(NY+1)
        X(i1,:) = X(i1-1,:);
    end

    for i2 = 2:(NX+1)
        Y(:,i2) = Y(:,i2-1);
    end

    grid on
    axis equal
    plot(X,Y,'r*')
end

% SOLVE THE PDE

% CONSTRUCT THE 5 DIAGONALS

% D1

D1 = zeros(JMAX*IMAX,1);
D1 = D1 + (-2*(DX^2+DY^2));

for i1 = 1:JMAX
    D1(i1) = 1; D1(i1+(IMAX-1)*JMAX) = 1; % exclude boundary nodes on the right and left line
end

for i1 = (JMAX+1):JMAX:((IMAX-2)*JMAX+1)
    D1(i1) = 1;  D1(i1+JMAX-1) = 1; % exclude boundary nodes on the up and down horizontal lines (above and below internal nodes) 
end

% D2 & D3

for i1 = 2:(JMAX*IMAX)
    
    if D1(i1) < 0
        D2(i1-1) = DX^2; 
    else
        D2(i1-1) = 0;
    end

    if D1(i1-1) < 0
        D3(i1-1) = DX^2; 
    else
        D3(i1-1) = 0; 
    end
end

% D4 & D5

for i1 = (JMAX+1):(IMAX*JMAX)

    if D1(i1) < 0
         D4(i1-JMAX) = DY^2;
    else
         D4(i1-JMAX) = 0;
    end
    if D1(i1-JMAX) < 0
        D5(i1-JMAX) = DY^2;
    else
        D5(i1-JMAX) = 0;
    end
end

% B VECTOR

for i1 = 1:(IMAX*JMAX)
    B(i1) = 4*DX^2*DY^2;
end

for i1 = 1:JMAX
    B(i1) = (((i1-1)*DY)^2); % LEFT BOUNDARY
    B(i1+(IMAX-1)*JMAX) = (((i1-1)*DY)^2) + (NX*DX)^2; % RIGHT BOUNDARY
end

n = 0;

for i1 = (JMAX+1):JMAX:((IMAX-2)*JMAX+1)
    n = n + 1;
    B(i1) = (n*DX)^2; %BOTTOM(BOUNDARY
    B(i1+JMAX-1) = (n*DX)^2 + (NY*DY)^2; %TOP BOUNDARY
end

%B = B';

% JACOBI Solution

w = input('Give the relaxation term (w): '); % set relaxation

X = B;
X_star = X;
X_new = X;

Ro_2 = 0;
R_2 = 0;
R_Ro = 10000;

for i1 = 1:(IMAX*JMAX)
    if D1(i1) < 0 
        n = n+1;
        X_new(i1) = 0.85; % random value for first iteration (only for internal nodes)
        %Ro_2 = (B(i1)-X_new(i1)*D1(i1)-X_new(i1-JMAX)*D4(i1-JMAX)-X_new(i1-1)*D2(i1-1)-X_new(i1+1)*D3(i1)-X_new(i1+JMAX)*D5(i1))^2 + Ro_2;
    end
end

for i1 = 1:(IMAX*JMAX)
    if D1(i1) < 0
        Ro_2 = (B(i1)-X_new(i1)*D1(i1)-X_new(i1-JMAX)*D4(i1-JMAX)-X_new(i1-1)*D2(i1-1)-X_new(i1+1)*D3(i1)-X_new(i1+JMAX)*D5(i1))^2 + Ro_2;
    end
end

X = X_new;
Ro = sqrt(Ro_2);

% ITERATIONS

figure;
h = animatedline('Color', 'b', 'LineWidth', 2);
addpoints(h, 0, 1);

set(gca, 'YScale', 'log');
xlabel('Iteration');
ylabel('Residual (R/Ro)');
grid on;

n = 0;

while R_Ro > 10^(-6)

for i = 1:(JMAX*IMAX)
    if D1(i) < 0
        X_star(i) = 1/D1(i)*(B(i)-X(i-JMAX)*D4(i-JMAX)-X(i-1)*D2(i-1)-X(i+1)*D3(i)-X(i+JMAX)*D5(i));
        X_new(i) = w*X_star(i) + (1-w)*X(i);
    end
end

for i = 1:(IMAX*JMAX)
    if D1(i) < 0
       R_2 = (B(i)-X_new(i)*D1(i)-X_new(i-JMAX)*D4(i-JMAX)-X_new(i-1)*D2(i-1)-X_new(i+1)*D3(i)-X_new(i+JMAX)*D5(i))^2 + R_2;
    end
end

n = n + 1;
R_Ro = sqrt(R_2)/Ro;
X = X_new;
R_2 = 0;
ITER = n;

addpoints(h, ITER, R_Ro);
drawnow limitrate;
end

% PLOT ISOLINES (CONTOUR)

x = 0:DX:1;
y = 0:DY:1;

[XX, YY] = meshgrid(x, y);

U = reshape(X, JMAX, IMAX);

figure;
contour(XX, YY, U, 30, 'ShowText', 'on');

colorbar;
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('Isolines of Numerical Solution');