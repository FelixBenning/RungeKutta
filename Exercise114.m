%% Test RungeKutta

A=[5/12,-1/12; 3/4,1/4];
c=[1/3;1];
b=[3/4,1/4]; 

rk=RungeKutta(A,b,c);


odefun = @(t,y) [3/4*y(1)-1/2*y(1)*y(2); -y(2)+y(1)*y(2)];

[T,Y,iterations]=rk.implWithFixedPointIter(odefun, 0:1/10:10, [1;1]);
plot(Y(1,:),Y(2,:))
xlabel("x")
ylabel("y")
%% iterations comparison
odefun = @(t,y) [3/4*y(1)-1/2*y(1)*y(2); -y(2)+y(1)*y(2)];

implEuler=RungeKutta(1,1,1);
[TEuler,~,EulerIter]=implEuler.implWithFixedPointIter(odefun, 0:1/10:10, [1;1]);
[TEulNew, ~, EulNewIter]=implEuler.implWithNewton(odefun, 0:1/10:10, [1;1]);
[TEulNewSimpl, ~, EulNewIterSimpl]=implEuler.implWithSimplNewton(odefun, 0:1/10:10, [1;1]);

implMidpoint=RungeKutta(1/2,1,1/2);
[TMidpoint,~,MidIter]=implMidpoint.implWithFixedPointIter(odefun, 0:1/10:10, [1;1]);
[TMidNew, ~, MidNewIter]=implMidpoint.implWithNewton(odefun, 0:1/10:10, [1;1]);
[TMidNewSimpl, ~, MidNewIterSimpl]=implMidpoint.implWithSimplNewton(odefun, 0:1/10:10, [1;1]);

plot(TEuler,EulerIter, 'o')
hold on
plot(TEulNew,EulNewIter, 'o')
plot(TEulNewSimpl,EulNewIterSimpl, 'o')
plot(TMidpoint,MidIter, 'o')
plot(TMidNew,MidNewIter, 'o')
plot(TMidNewSimpl,MidNewIterSimpl, 'o')
hold off
xlabel("t")
ylabel("iterations")
legend("Euler fixedPoint", "Euler Newton", "Euler Simplified Newton","Midpoint fixedPoint", "Midpoint Newton", "Midpoint Simplified Newton")