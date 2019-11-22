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
implEuler=RungeKutta(1,1,1);
implMidpoint=RungeKutta(1/2,1,1/2);

odefun = @(t,y) [3/4*y(1)-1/2*y(1)*y(2); -y(2)+y(1)*y(2)];
[TEuler,~,EulerIter]=implEuler.implWithFixedPointIter(odefun, 0:1/10:10, [1;1]);
[TMidpoint,~,MidIter]=implMidpoint.implWithFixedPointIter(odefun, 0:1/10:10, [1;1]);
plot(TEuler,EulerIter, 'o')
hold on
plot(TMidpoint,MidIter, 'o')
hold off
xlabel("t")
ylabel("iterations")
legend("Euler","Midpoint")