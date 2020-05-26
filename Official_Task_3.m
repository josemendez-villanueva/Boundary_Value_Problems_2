N = 800;
M = 800;
a = 1;
d = .1;



xx = linspace(0, 1, N+1)';
Insidex = xx(1:end-1);
deltax = 1/N;

Pe = abs(a/d)*deltax;

tf = 1;
deltat = tf/M;
tt = linspace(0, tf, M+1);
[X, T] = meshgrid(xx, tt);

Solution = zeros(N, M+1);

ICeq = exp(-100*(Insidex - 0.5).^2);

Solution(:, 1) = ICeq;
uold = ICeq;


for i = 1:M
    unew = ConvectionDiffusionEq(uold, a, d, deltat);
    Solution(:, i + 1) = unew;
    uold = unew;
end

Solution = [Solution; Solution(1, :)];

figure(1);
surf(X, T, Solution', 'Edgecolor', 'none');
title('Convection-Diffusion Eq Pe=.1');
xlabel('X');
ylabel('T');

display(Pe)

function unew = ConvectionDiffusionEq(u, a, d, dt)
    N = length(u);
    M = 1/dt;
    deltax = 1/N;
    deltaxSQRD = deltax^2;
    
    Sub = d/deltaxSQRD + a/(2*deltax);
    Main = - 2*d/deltaxSQRD;
    Sup = d/deltaxSQRD - a/(2*deltax);
    A = diag(Sub*ones(N-1,1),-1) + diag(Main*ones(N,1) ,0) + diag(Sup*ones(N-1,1) ,1);
    A(1,N) = Sub;
    A(N,1) = Sup;
    unew = TRAP(A, u, dt);
end