function [Ga] = get_g0_single(a,type_norm)

    % Function that takes an auto correlation and returns the Normalization
    % Value G0 (0 delay), default is to use the G0_intp, interpolate from 
    % G1, G2, G3. Other options are to use the G0 (would include shot noise), 
    % or use G1.

    switch type_norm
        case 'none'
            Ga=1;
    
        case 'G1'
            Ga=a(2);
    
        case 'G0'
            Ga=a(1);
    
        case 'G0_intp'
            Ga=a(2) - (a(3)-a(2));
    
    end
end