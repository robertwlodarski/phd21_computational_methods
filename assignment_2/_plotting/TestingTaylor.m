% Testing Taylor approximations

f0      = @(x) exp(x);
f1      = @(x) 1+x;
f2      = @(x) 1+x+x^2/2;
f3      = @(x) 1+x+x^2/2+x^3/6;

figure(1);
fplot(f0,[0,1],'LineWidth',2.5);
hold on;
fplot(f1,[0,1],'LineWidth',2.5);
fplot(f2,[0,1],'LineWidth',2.5);
fplot(f3,[0,1],'LineWidth',2.5);
grid on;
legend('$\exp(x)$','$1+x$','$1+x+\frac{x^2}{2}$','$1+x+\frac{x^2}{2}+\frac{x^3}{6}$', 'Location', 'best','interpreter','latex');
hold off;


% Testing 3rnd degree Taylor for Gumbel trick 
A       = 10;
B       = 13;

g0      = @(x) exp(A/x) / (exp(A/x)+exp(B/x));
g1      = @(x) 1 / (1 + exp(B-A)) + (B-A) / ((1 + exp(B-A))^2) * (x-1);

g0(1) 
g1(1)

g0(0.1)
g1(0.1)

g0(0.01)
g1(0.01)