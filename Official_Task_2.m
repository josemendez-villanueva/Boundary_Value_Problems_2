N = 800;
M = 800;

xx = linspace(0, 1, N+1)';
Insidex = xx(1:end-1);
deltax = 1/N;
Courant = .9;

tf = 5;
deltat = tf/M;
tt = linspace(0, tf, M+1);
[X, T] = meshgrid(xx, tt);

Solutions = zeros(N, M+1);
ICeq = exp(-100*(Insidex - 0.5).^2);
Solutions(:, 1) = ICeq;
uOld = ICeq;

A = zeros(N, 1);
A(1) = -2; A(2)=1;
BigA = toeplitz(A)/deltax^2;

NORM = zeros(M+1, 1);
NORM(1) = rms(ICeq);

for i = 1:M
    unew = LAXWENDROF(uOld, Courant);
    Solutions(:, i + 1) = unew;
    NORM(i + 1) = rms(unew);
    uOld = unew;
end

Solutions = [Solutions; Solutions(1, :)];

figure(1);
surf(X, T, Solutions','EdgeColor','none');
title('Advection equation');
xlabel('X');
ylabel('T');

figure(2);
plot(tt, NORM);
title('Norm of Solution');
xlabel('T');
ylabel('Norm');

function unew = LAXWENDROF(u, Courant)
    Lower = (Courant/2)*(1+Courant);
    Main = 1 - Courant^2;
    Upper = -(Courant/2)*(1-Courant);
    N = length(u);
    MainMatrix = diag(Lower*ones(N-1,1),-1) + diag(Main*ones(N,1),0) + diag(Upper*ones(N-1,1),1);
    MainMatrix(N,1) = Upper;
    MainMatrix(1,N) = Lower;
    
    unew = MainMatrix*u;
end
