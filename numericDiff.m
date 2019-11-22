function[df]=numericDiff(func,x,TOL)
    nn=length(x);
    pdiff = @(i) directionalDerivative(func, x, i, TOL);
    id=eye(length(x));
    df=splitapply(pdiff,id,1:nn);
end
function[dfDir]=directionalDerivative(func,x,direction,TOL)
    eps=10^-6;
    delta=1;
    dfDirNext=(func(x+eps*direction)-func(x))/eps;
    while delta>TOL
        eps=eps/2;
        dfDirPrev=dfDirNext;
        dfDirNext=(func(x+eps*direction)-func(x))/eps;
        delta=norm(dfDirPrev-dfDirNext,inf);
    end
    dfDir=dfDirNext;
end

