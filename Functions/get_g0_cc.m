function [Ga,Gb,Gc] = get_g0_cc(a,b,c,type_norm)

switch type_norm
    case 'none'
        Ga=1;
        Gb=1;
        Gc=1;
    case 'Max'
        Ga=max(a);
        Gb=max(b);
        Gc=max(c);
    case 'G0'
        Ga=a(ceil(length(a)/2));
        Gb=b(ceil(length(b)/2));
        Gc=c(ceil(length(c)/2));

end
end