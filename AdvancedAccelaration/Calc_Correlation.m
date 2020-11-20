function [step_list,ans_list] = Calc_Correlation(sample)
    sample = sample - mean(sample);
    sample = sample/sqrt(var(sample));
    N = length(sample);
    N_output = 100;
    step = fix(N/(2*N_output));
    step_list = step:step:N_output*step;
    ans_list = zeros([N_output,1]);
    for step_index = 0:1:(N_output-1)
        distance = step_index * step;
        temp = sample(1:(N-distance)).*sample((1+distance):N);
        ans_list(step_index+1) = mean(temp);
    end

end

