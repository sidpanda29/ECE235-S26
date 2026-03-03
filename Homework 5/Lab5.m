% Sidorian Pandovski
% ECE 235 
% Homework 5
%
% a) - Homework 4
% b) - Homework 4
% c) - Homework 4
%
% d) On each graph from c), in a color of your choice, shade the area under
% the probability-density curve that corresponds to the range of positions 
% between 〈𝑥〉 − ∆𝑥 and 〈𝑥〉 + ∆𝑥 . (Hint: Use the MATLAB area command or 
% the fill command;).
% 
% e) Calculate the probability of finding the particle in the range between
%〈𝑥〉 − ∆𝑥 and 〈𝑥〉 + ∆𝑥. This number is equal to the area of the pretty 
% shaded region from part d. For wave functions 1 and 2, do this 
% analytically. For wave function 3, do it using MATLAB’s integral command 
% or piecewise rectangular integration (as we did in HW 2).
%
% f) - Homework 4
%
% g) Calculate the product of the uncertainties ∆𝑥∆𝑝. Does it satisfy the 
% Heisenberg uncertainty relation ∆𝑥∆𝑝 ≥ ℏ/2 ?
%
% test wave functions:
% 1) Psi(x,0) = Ax(x-L) where L is known. the particle is confined to the
% region [0,L] -> therefore zero outside this region
%
% 2) Psi(x,0) = Asin(n*pi*x/L) where n is an arbitrary positive integer and
% L is known. the particle is confined to the region [0,L] -> therefore
% zero outside this region
%
% 3) Psi(x,0) = Aexp[(-(x-xNaught)^2)/sigma^2 + i*kNaught*x] the particle
% can move in all 1D space

% Defining variables:

eMass   = 9.1E-31; % kilograms
L       = 10E-9; % meters
sigma   = 1E-10; % meters
x0      = 2E-10; % meters
k0      = 2E9;% radians/meter
n       = [1 2 5]; 

x = linspace(0, L, 1000);
x3 = linspace(x0-2*sigma, x0+2*sigma, 1000);

% Normalization values:

A1 = sqrt(30/L^5);
A2 = sqrt(2 / L);
A3 = (2 / (pi * sigma^2))^(1/4);

% Given wavefunctions:

Psi1 = A1 .* x .* (x-L);
Psi2 = A2 .* sin(((pi.*n')/L).*x);
Psi3 = A3 .* exp(((-(x3-x0).^2)/(sigma^2)) + (j * k0 .* x3));

figure;
plot(x, abs(Psi1).^2);
xlabel('x (m)'); ylabel('|\Psi|^2'); title('PDF: Wave Function 1');

figure;
plot(x, abs(Psi2).^2);
xlabel('x (m)'); ylabel('|\Psi|^2'); title('PDF: Wave Function 2 (n = 1,2,5)');

figure;
plot(x3, abs(Psi3).^2);
xlabel('x (m)'); ylabel('|\Psi|^2'); title('PDF: Wave Function 3');

%% Part d) Shade area under PDF between <x>-Dx and <x>+Dx

% --- Wave Function 1 ---
% <x> = L/2 by symmetry; <x^2> = 2L^2/7; (Dx)^2 = 2L^2/7 - L^2/4 = L^2/28
xMean1 = L / 2;
dx1    = L / sqrt(28);
xlo1   = xMean1 - dx1;
xhi1   = xMean1 + dx1;

mask1  = (x >= xlo1) & (x <= xhi1);
xFill1 = x(mask1);
yFill1 = abs(A1 .* xFill1 .* (xFill1 - L)).^2;

figure(1); hold on;
fill([xFill1(1), xFill1, xFill1(end)], [0, yFill1, 0], 'c', ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold off;

% --- Wave Function 2 (n = 1, 2, 5) ---
% <x> = L/2 by symmetry; <x^2> = L^2/3 - L^2/(2n^2*pi^2)
% (Dx)^2 = L^2/12 - L^2/(2n^2*pi^2)
xMean2    = L / 2;
shadeCols = {'m', 'g', 'r'};

figure(2); hold on;
for k = 1:length(n)
    nk   = n(k);
    dx2  = L * sqrt(1/12 - 1/(2 * nk^2 * pi^2));
    xlo2 = xMean2 - dx2;
    xhi2 = xMean2 + dx2;

    mask2  = (x >= xlo2) & (x <= xhi2);
    xFill2 = x(mask2);
    yFill2 = abs(A2 .* sin(nk * pi .* xFill2 / L)).^2;
    fill([xFill2(1), xFill2, xFill2(end)], [0, yFill2, 0], shadeCols{k}, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
hold off;

% --- Wave Function 3 ---
% |Psi3|^2 = A3^2 * exp(-2*(x-x0)^2/sigma^2)  [Gaussian centered at x0]
% <x> = x0;  (Dx)^2 = sigma^2/4  =>  Dx = sigma/2
xMean3 = x0;
dx3    = sigma / 2;
xlo3   = xMean3 - dx3;
xhi3   = xMean3 + dx3;

mask3  = (x3 >= xlo3) & (x3 <= xhi3);
xFill3 = x3(mask3);
yFill3 = A3^2 .* exp(-2 * (xFill3 - x0).^2 / sigma^2);

figure(3); hold on;
fill([xFill3(1), xFill3, xFill3(end)], [0, yFill3, 0], 'g', ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold off;

%% Part e) Probability of finding particle in [<x>-Dx, <x>+Dx]

% --- Wave Function 1 (analytical) ---
% int |Psi1|^2 dx = A1^2 * int x^2*(x-L)^2 dx
%                = A1^2 * int (x^4 - 2L*x^3 + L^2*x^2) dx
% Antiderivative: A1^2 * (x^5/5 - L*x^4/2 + L^2*x^3/3)
F1 = @(t) A1^2 .* (t.^5/5 - L.*t.^4/2 + L^2.*t.^3/3);
P1 = F1(xhi1) - F1(xlo1);

fprintf('=== Wave Function 1 ===\n');
fprintf('  <x> = %.4e m,  Dx = %.4e m\n', xMean1, dx1);
fprintf('  P(<x> +/- Dx) = %.4f\n\n', P1);

% --- Wave Function 2 (analytical, each n) ---
% int |Psi2|^2 dx = A2^2 * int sin^2(n*pi*x/L) dx
%                = A2^2 * (x/2 - L*sin(2*n*pi*x/L) / (4*n*pi))
fprintf('=== Wave Function 2 ===\n');
for k = 1:length(n)
    nk   = n(k);
    dx2  = L * sqrt(1/12 - 1/(2 * nk^2 * pi^2));
    xlo2 = xMean2 - dx2;
    xhi2 = xMean2 + dx2;

    F2 = @(t) A2^2 .* (t/2 - L .* sin(2*nk*pi.*t/L) / (4*nk*pi));
    P2 = F2(xhi2) - F2(xlo2);
    fprintf('  n=%d:  <x>=%.4e m,  Dx=%.4e m,  P=%.4f\n', nk, xMean2, dx2, P2);
end
fprintf('\n');

% --- Wave Function 3 (numerical via integral) ---
pdf3 = @(t) A3^2 .* exp(-2*(t - x0).^2 / sigma^2);
P3   = integral(pdf3, xlo3, xhi3);

fprintf('=== Wave Function 3 ===\n');
fprintf('  <x> = %.4e m,  Dx = %.4e m\n', xMean3, dx3);
fprintf('  P(<x> +/- Dx) = %.4f\n\n', P3);

