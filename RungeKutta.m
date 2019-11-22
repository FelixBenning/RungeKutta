classdef RungeKutta    
    properties
        %Butcher Array:
        A(:,:) {mustBeReal, mustBeFinite}
        b(1,:) {mustBeReal, mustBeFinite}
        c(:,1) {mustBeReal, mustBeFinite}
    end
    
    methods
        function obj = RungeKutta(A,b,c)
            obj=obj.setButcher(A,b,c);
        end
        
        function obj=setButcher(obj,A,b,c)
            ASize=size(A);
            if(ASize(1)~=length(c))
                throw(MException("A and c don't have matching dimensions"))
            end
            if(ASize(2)~=length(b))
                throw(MException("A and b don't have matching dimensions"))
            end
            obj.A=A;
            obj.b=b;
            obj.c=c;
        end
        
        function [T,Y, iterations]=implWithFixedPointIter(obj, odefun, tspan, y0)
            %Solve ode using the Fixed Point Iteration method to solve the
            %fixed point equation of the implicit RK scheme
            solver=@(func, startPoint) fixedPointIter(func, startPoint, 10^(-10), 10^6);
            [T,Y, iterations] = obj.implWithCustomSolver(solver,odefun,tspan,y0);
        end
        
        function [T,Y, iterations]=implWithNewton(obj, odefun, tspan, y0)
            %Solve ode using the Newton method to solve the
            %fixed point equation of the implicit RK scheme
            function [fixedPoint, hasConverged, iterations]=solver(func, startPoint)
                zeroProblem = @(var) var-func(var);
                [fixedPoint, hasConverged, iterations]=newtonSolver(zeroProblem, startPoint, 10^(-10), 10^6);
            end
            [T,Y, iterations] = obj.implWithCustomSolver(@solver,odefun,tspan,y0);
        end
        
        function [T,Y, iterations]=implWithSimplNewton(obj, odefun, tspan, y0)
            %Solve ode using the Simplified Newton method to solve the
            %fixed point equation of the implicit RK scheme
            function [fixedPoint, hasConverged, iterations]=solver(func, startPoint)
                zeroProblem = @(var) var-func(var);
                [fixedPoint, hasConverged, iterations]=simplNewtonSolver(zeroProblem, startPoint, 10^(-10), 10^6);
            end
            [T,Y, iterations] = obj.implWithCustomSolver(@solver,odefun,tspan,y0);
        end
        
        function [T,Y, iterations]=implWithCustomSolver(obj, solver, odefun, tspan, y0)
            %solver is called to find the fixed point of the implicit RK equation
            %solver needs to be a function with the signature
            %[fixedPoint, hasConverged, iterations] = solver(function, startPoint) 
            s=length(obj.c);
            
            if length(tspan)==2
                tt = linspace(tspan(1), tspan(2), 1000);
            else
                tt = tspan;
            end

            yy = zeros(length(y0), length(tt));
            iterations = zeros(1, length(tt));
            z0=zeros(length(y0), s);

            yy(:,1)= y0;
            for ii = 1:(length(tt)-1)
               tau=(tt(ii+1)-tt(ii));
               if s>1
                   %splitapply(odefun,t,z,(1:s)) returns [odefun(t(i),z(i))]_{i=1:s} 
                   %for entries t(i) of t, and columns z(i) of z 
                   func=@(var) tau*obj.A*(splitapply(odefun, tt(ii)+tau*obj.c', yy(:,ii)+var, 1:s))';
               else
                   func=@(var) tau*obj.A*odefun(tt(ii)+tau*obj.c, yy(:,ii)+var);
               end
               %%%%%%% main step %%%%%%%%%%%
               [zz, hasConverged, iterations(1,ii)]=solver(func, z0);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if ~hasConverged
                   throw(MException("Solver could not find fixed point"))
               end
               if s>1
                   yy(:,ii+1)=yy(:,ii)+tau*(splitapply(odefun, tt(ii)+tau*obj.c', yy(:,ii)+zz, (1:s)))*obj.b';
               else
                   yy(:,ii+1)=yy(:,ii)+tau*obj.A*odefun(tt(ii)+tau*obj.c, yy(:,ii)+zz);
               end
            end

            T=tt;
            Y=yy;
        end
    end
end

