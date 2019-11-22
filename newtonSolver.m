function [fixedPoint, hasConverged, cycles] = newtonSolver(func, startPoint, TOL, patience)
    df=numericDiff(func,startPoint, 0.5*TOL);
    prev=startPoint-df\func(startPoint);
    prevDelta=norm(prev-startPoint, inf);
    
    cycles=1;
    
    while 1
        df=numericDiff(func,prev, 0.5*TOL);
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