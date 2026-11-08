% Finite Difference Method for 1 x 1 m^2 for POISSON equation ( = 4 )

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




% B VECTOR & 5 DIAGONALS

% INTERNAL NODES

for i = 1:IMAX
    for j = 1:JMAX
        k = JMAX*(i-1) + j;
        D1(k) = -2*(DX^2 + DY^2);
        D2(k) = DX^2;
        D3(k) = DX^2;
        D4(k) = DY^2;
        D5(k) = DY^2;
        B(k) = 4*DX^2*DY^2;
    end
end

%LEFT BOUNDARY

for i = 1:IMAX:IMAX
    for j = 1:JMAX
        k = JMAX*(i-1) + j;
        D1(k) = 1;
        D2(k) = 0;
        D3(k) = 0;
        D4(k) = 0;
        D5(k) = 0;
        B(k) = (((j-1)*DY)^2);
    end
end

%RIGHT BOUNDARY

for i2 = 1:IMAX:IMAX
    for j = 1:JMAX
        i = IMAX - i2 + 1;
        k = JMAX*(i-1) + j;
        D1(k) = 1;
        D2(k) = 0;
        D3(k) = 0;
        D4(k) = 0;
        D5(k) = 0;
        B(k) = ((j-1)*DY)^2 + (NX*DX)^2;
    end
end

%BOTTOM BOUNDARY

for i = 1:IMAX
    for j = 1:JMAX:JMAX
        k = JMAX*(i-1) + j;
        D1(k) = 1;
        D2(k) = 0;
        D3(k) = 0;
        D4(k) = 0;
        D5(k) = 0;
        B(k) = ((i-1)*DX)^2;
    end
end

%TOP BOUNDARY

for i = 1:IMAX
    for j2 = 1:JMAX:JMAX
        j = JMAX - j2 +1;
        k = JMAX*(i-1) + j;
        D1(k) = 1;
        D2(k) = 0;
        D3(k) = 0;
        D4(k) = 0;
        D5(k) = 0;
        B(k) = ((i-1)*DX)^2 + ((JMAX-1)*DY)^2;
    end
end



% SOLVE THE PDE


% JACOBI Solution

w = input('Give the relaxation term (w): '); % set relaxation

X = B;
X_star = X;
X_new = X;

Ro_2 = 0; %for measuring the square terms ofn1st residual
R_2 = 0; %for measuring the square terms of residual
R_Ro = 10000; %first big normalized residual for while check to keep running before producing and R_Ro


for i = 2:(IMAX-1)
    for j = 2:(JMAX-1)
        k = JMAX*(i-1) + j;
        X_new(k) = 0.85; % random value for first iteration (only for internal nodes)
    end
end

for i = 2:(IMAX-1)
    for j = 2:(JMAX-1)
        k = JMAX*(i-1) + j;
        Ro_2 = (B(k)-X_new(k)*D1(k)-X_new(k-JMAX)*D4(k)-X_new(k-1)*D2(k)-X_new(k+1)*D3(k)-X_new(k+JMAX)*D5(k))^2 + Ro_2;;
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
max_iter = 1000000;

while R_Ro > 10^(-6) && n < max_iter % Drop the initial residual by 6 orders of magnitude or until max number of iterations


    for i = 2:(IMAX-1)
        for j = 2:(JMAX-1)
            k = JMAX*(i-1) + j;
            X_star(k) = 1/D1(k)*(B(k)-X(k-JMAX)*D4(k)-X(k-1)*D2(k)-X(k+1)*D3(k)-X(k+JMAX)*D5(k));
            X_new(k) = w*X_star(k) + (1-w)*X(k);
        end
    end

    for i = 2:(IMAX-1)
        for j = 2:(JMAX-1)
            k = JMAX*(i-1) + j;
            R_2 = (B(k)-X_new(k)*D1(k)-X_new(k-JMAX)*D4(k)-X_new(k-1)*D2(k)-X_new(k+1)*D3(k)-X_new(k+JMAX)*D5(k))^2 + R_2;
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

x_grid = 0:DX:1;
y_grid = 0:DY:1;

[XX, YY] = meshgrid(x_grid, y_grid);

U = reshape(X, JMAX, IMAX);

figure;
contour(XX, YY, U, 30, 'ShowText', 'on');

colorbar;
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('Isolines of Numerical Solution');