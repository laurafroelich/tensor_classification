function [sigma, corestd, coreNoiseStd, UNoiseStd] = get_config(signal_level)

switch signal_level
    case 'low'
        sigma = 4;
        corestd = 1;
        coreNoiseStd = sqrt(3);
        UNoiseStd = sqrt(0.01);
        
    case 'medium'
        sigma = 5;
        corestd = sqrt(3);
        coreNoiseStd = sqrt(3);
        UNoiseStd = sqrt(0.01);
        
    case 'high'
        sigma = 5;
        corestd = sqrt(5);
        coreNoiseStd = 1;
        UNoiseStd = sqrt(0.01);
        
    otherwise
        disp('no such level')

end