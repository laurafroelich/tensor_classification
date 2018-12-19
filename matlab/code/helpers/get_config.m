function [sigma, corestd, coreNoiseStd] = get_config(signal_level)

switch signal_level
    case 'low'
        sigma = 25;
        corestd = 0.25;
        coreNoiseStd = sqrt(3);
        
    case 'medium'
        sigma = 25;
        corestd = 0.5;
        coreNoiseStd = sqrt(3);
        
    case 'high'
        sigma = 25;
        corestd = 1;
        coreNoiseStd = 1;
        
    otherwise
        disp('no such level')

end