function tests = test_gradients
    tests = functiontests(localfunctions);
end

function test_parafac_lda_gradient(testCase)
    analytical_vs_numerical_gradient(testCase, 'U1',...
        @parafacldaobj_matrixdata)
    analytical_vs_numerical_gradient(testCase, 'U2',...
        @parafacldaobj_matrixdata)
end

function test_tucker_lda_gradient(testCase)
    analytical_vs_numerical_gradient(testCase, 'U1',...
        @tensorsldaobj_matrixdata)
    analytical_vs_numerical_gradient(testCase, 'U2',...
        @tensorsldaobj_matrixdata)
end

function analytical_vs_numerical_gradient(testCase, matrix_name, ...
    objective_function)
    p = 5;
    q = 7;
    nsamples = 100;
    [Xs, ys] = get_data(nsamples, [p, q]);
    
    [classmeandiffs, observationdiffs, nis] = ...
        classbased_differences(Xs, ys);
    
    k1 = 3;
    k2 = 3;
    myfun = @(U) objective_function(U,...
        classmeandiffs, observationdiffs, nis, k2, k2);

    initial_matrices_struct = get_initial_matrices(p, q, k1, k2);
    numerically_estimated_gradient = get_numerically_estimated_gradient(...
        myfun, initial_matrices_struct, matrix_name);
    
    [F, analyticalderiv] = myfun(initial_matrices_struct);
    
    actual_difference = norm(...
        numerically_estimated_gradient-...
        analyticalderiv.(matrix_name), 'fro')/...
        norm(analyticalderiv.(matrix_name), 'fro');
    
    verifyEqual(testCase, 0, actual_difference, 'AbsTol', 1e-5)
end

function matrices_in_struct = get_initial_matrices(p, q, k1, k2)

    if k1 ~= k2
        warning('k2 must be same value as k1, will be set as such')
    end
    
    [U, ~, ~] = svd(randn(p));
    U1 = U(:, 1:k1);

    k2 = k1;
    [U, ~, ~] = svd(randn(q));
    U2 = U(:, 1:k2);

    matrices_in_struct.U1 = U1;
    matrices_in_struct.U2 = U2;
end

function [Xs, classes] = get_data(nsamples, xdims)

    if nargin < 2
        xdims = [5, 7];

        if nargin < 1
            nsamples = 500;
        end
    end

    classes = round(rand(nsamples,1))+1;
    Xs = cell(1, nsamples);

    for isample = 1:nsamples
        if classes(isample) == 1
            Xs{isample} =  randn(xdims);
        else

            Xs{isample} =  randn(xdims) + 10;
        end
    end

    
end

function estderiv = get_numerically_estimated_gradient(myfun, ...
    Ustruct, matrix_name)

    [I, J] = size(Ustruct.(matrix_name));
    myeps = 1e-9;
    estderiv = zeros(I, J);

    for i=1:I
        for iparam=1:J
            current_u = Ustruct.(matrix_name);
            wp = current_u;
            wm = current_u;
            wp(i, iparam) = wp(i, iparam) + myeps;
            wm(i, iparam) = wm(i, iparam) - myeps;
            Ustruct.(matrix_name) = wp;

            Fp = myfun(Ustruct);

            Ustruct.(matrix_name) = wm;
            Fm = myfun(Ustruct);
            estderiv(i, iparam) = (Fp-Fm)/(2*myeps);
        end
    end

end

