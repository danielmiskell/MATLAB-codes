%  D. Hewett
%  MATH3603 Numerical Methods
%  Demo 1 - nonlinear equations
%% Setup
close all, clear all, clc
tol=1e-8;
nmax=100;
fs=16;
set(0,'defaulttextfontsize',fs);
set(0,'defaultaxesfontsize',fs);
%% Bisection
format long, format compact
f=@(x)sin(2*x)-1+x;
x=linspace(-3,3,101);
plot(x,f(x),'LineWidth',2)
grid on
xlabel('x')
ylabel('f(x)')
hold on
[zero,res,niter,itersb]=bisection(f,-1,1,tol,nmax)
%scatter(zero,res)
for i=1:10
    scatter(itersb(i,1),f(itersb(i,1)),'r')
    pause
end
% Play with nmax, tol, a, b
%% Fixed point iterations
phi1=@(x)1-sin(2*x);
phi2=@(x).5*asin(1-x);
figure
x=linspace(0,2,101);
plot(x,phi1(x),'-k','LineWidth',2)
hold on
plot(x,phi2(x),'--b','LineWidth',2)
plot(x,x,'-.r','LineWidth',2)
grid on
xlabel('x')
legend('\phi_1(x)','\phi_2(x)','x','Location','NorthWest')
[fixp,res,niter,itersfp1]=fixpoint(phi1,0.7,tol,nmax)
%scatter(fixp,fixp,'x')
[fixp,res,niter,itersfp2]=fixpoint(phi2,0.7,tol,nmax)
%scatter(fixp,fixp)
for i=1:niter+1
    scatter(itersfp2(i,1),itersfp2(i,1),'r')
    pause
end
%% Newton method
df = @(x)2*cos(2*x)+1;
[zero,res,niter,itersn]=newton(f,df,0.7,tol,nmax)
%% Chord method
[zero,res,niter,itersc]=chord(f,-1,1,0.7,tol,nmax)
%% Secant method
x1=chord(f,-1,1,0.7,tol,1)
[zero,res,niter,iterss]=secant(f,0.7,x1,tol,nmax)
%% Comparison of convergence
zeroex=newton(f,df,0.7,1e-15,nmax)
errb=abs(itersb-zeroex);
errfp2=abs(itersfp2-zeroex);
errn=abs(itersn-zeroex);
errc=abs(itersc-zeroex);
errs=abs(iterss-zeroex);
figure
semilogy(errb,'-k','LineWidth',2)
hold on
semilogy(errfp2,'--b','LineWidth',2)
semilogy(errn,':m','LineWidth',2)
semilogy(errc,'-.r','LineWidth',2)
semilogy(errs,'-g','LineWidth',2)
xlabel('iteration')
ylabel('error')
legend('bisection','fixed point \phi_2','Newton','chord','secant')
%% Multi-dim Newton
f1 = @(x1,x2)x1.^2-2*x1.*x2-2;
f2 = @(x1,x2)x1+x2.^2+1;
[x1,x2] = meshgrid(-6:.1:2,-4:.1:4);
figure
contour(x1,x2,f1(x1,x2),[0,0],'k'); hold on;
contour(x1,x2,f2(x1,x2),[0,0],'b--')
legend('f_1(x_1,x_2)=0','f_2(x_1,x_2)=0','Location','South')
f=@(x)[f1(x(1),x(2));f2(x(1),x(2))];
Jf = @(x)[2*x(1)-2*x(2),-2*x(1);1,2*x(2)];
[zero,res,niter,itersn1]=newton_multidim(f,Jf,[-1;0],tol,nmax)
for i=1:niter+1
    scatter(itersn1(i,1),itersn1(i,2),'r')
    pause
end
[zero,res,niter,itersn2]=newton_multidim(f,Jf,[-5;-3],tol,nmax)
for i=1:niter+1
    scatter(itersn2(i,1),itersn2(i,2),'b')
    pause
end
% [zero,res,niter,itersn3]=newton_multidim(f,Jf,[5;3],tol,nmax)
% for i=1:niter+1
%     scatter(itersn3(i,1),itersn3(i,2),'g')
% end
