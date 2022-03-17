function [Ga] = get_g0_cc_single(a,type_norm)
    % Function that takes a cross correlation and returns the Normalization
    % Value G0, default is to use the center value (0 delay)

    switch type_norm
        case 'none'
            Ga=1;
        case 'max'
            Ga=max(a);
    
        case 'G0'
            Ga=a(ceil(length(a)/2));
    end
end