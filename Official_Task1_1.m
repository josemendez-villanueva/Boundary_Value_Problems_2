N = 20;
M = 100;
xx = linspace(0, 1, N+2)';
Insidex = xx(2:end-1);
deltax = 1/(N+1);
tF = .11;
deltat = tF/M;
tt = linspace(0, tF, M+1);
[X, T] = meshgrid(xx, tt);

A = zeros(N,1);
A(1) =-2; A(2) =1;
BigA = toeplitz(A)/deltax^2;
Solution = zeros(N, M+1);
ICeq = exp(-50.*(Insidex-.3).^2);

Solution(:, 1) = ICeq;
uOld = ICeq;

for i = 1:M
    uNew = Euler(BigA, uOld, deltat);
    Solution(:, i+1) = uNew;
    uOld = uNew;
end

Solution = [zeros(1, length(tt)); Solution; zeros(1, length(tt))];

figure(1);
surf(X, T, Solution', 'Edgecolor', 'none')
title('Stable Diffusion equation @ tF = .11');
xlabel('X');
ylabel('T');

CFL = deltat/deltat^2;
display(CFL);

function unew = Euler(Tdx, uold, dt)
    unew = uold + dt*Tdx*uold;
end