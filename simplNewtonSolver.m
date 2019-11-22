function [fixedPoint, hasConverged, cycles] = simplNewtonSolver(func, startPoint, TOL, patience)
    df=differentiate(func,startPoint);
    prev=startPoint-df\func(startPoint);
    prevDelta=norm(prev-startPoint, inf);
    
    cycles=1;
    
    while 1
        next=prev-df\func(prev);
        delta=norm(next-prev,inf);
        theta=min(delta/prevDelta, 0.9);
        
        if ((delta*theta/(1-theta))<0.5*TOL || cycles>patience)
            break;
        end
        cycles=cycles+1;
        prev=next;
        prevDelta=delta;
    end
    if(cycles<patience)
        fixedPoint=next;
        hasConverged=true;
    else
        fixedPoint=[];
        hasConverged=false;
    end
end
