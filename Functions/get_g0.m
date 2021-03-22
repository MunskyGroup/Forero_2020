function [Ga,Gb,Gc] = get_g0(a,b,c,type_norm)

switch type_norm
    case 'none'
        Ga=1;
        Gb=1;
        Gc=1;
    case 'G1'
        Ga=a(2);
        Gb=b(2);
        Gc=c(2);
    case 'G0'
        Ga=a(1);
        Gb=b(1);
        Gc=c(1);
    case 'G0_intp'
        Ga=a(2) - (a(3)-a(2));
        Gb=b(2) - (b(3)-b(2));
        Gc=c(2) - (c(3)-c(2));
end
end