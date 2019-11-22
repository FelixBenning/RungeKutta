function[df]=differentiate(func,x,TOL)
    nn=length(x);
    df=splitapply(partdiff,func,x',1:nn,TOL,1:nn);
end
function[dfi]=partdiff(func,x,i,TOL)
    eps=0.1;
    delta=1;
    ei=eye(length(x),i);
    dfinext=(func(x+eps*ei)-func(x))/eps;
    while delta>TOL
        eps=eps/2;
        dfiprev=dfinext;
        dfinext=(func(x+eps*ei)-func(x))/eps;
        delta=norm(dfiprev-dfinext,inf);
    end
    dfi=dfinext;
end

