function [x,y,d]=extendedEuclidean(a,b)
    if b==0
        x=1;
        y=0;
        d=a;
        return;
    end
    [x1,y1,d1]=extendedEuclidean(b,mod(vpi(a),vpi(b)));  

    x=y1;
    y=vpi(x1)-floor(vpi(a)/vpi(b))*vpi(y1);
    d=d1;

end 