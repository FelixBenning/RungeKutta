function [fixedPoint, hasConverged, cycles, theta] = fixedPointIter(func, startPoint, TOL, patience)
    % finds fixpoints of Contractions by repeatedly applying func, starting
    % with start point. Stops when TOL criteria met or patience runs out
    prev=func(startPoint);
    prevDelta=norm(prev-startPoint, inf);
    next=func(prev);
    delta=norm(next-prev, inf);
    theta=min(delta/prevDelta, 0.9);
    cycles=1;
    while((delta*theta/(1-theta))>0.5*TOL && cycles<patience)
        prev=next;
        prevDelta=delta;
        next=func(prev);
        delta=norm(next-prev, inf);
        theta=min(delta/prevDelta, 0.9);
        cycles=cycles+1;
    end
    if(cycles<patience)
        fixedPoint=next;
        hasConverged=true;
    else
        fixedPoint=[];
        hasConverged=false;
    end
end