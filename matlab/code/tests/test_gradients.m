function tests = test_gradients            
    tests = functiontests(localfunctions);
end

function test_parafac_lda_gradient(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/manopt/', 'IncludeSubfolders', true));
        
    analytical_vs_numerical_gradient(testCase, 'U1',...
        @parafacldaobj_matrixdata)
    analytical_vs_numerical_gradient(testCase, 'U2',...
        @parafacldaobj_matrixdata)
     %analytical_vs_numerical_gradient(testCase, 'U3',...
     %    @parafacldaobj_matrixdata)
end

function test_tucker_lda_gradient(testCase)    
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));

    analytical_vs_numerical_gradient(testCase, 'U1',...
        @tensorsldaobj_matrixdata)
    analytical_vs_numerical_gradient(testCase, 'U2',...
        @tensorsldaobj_matrixdata)
    %analytical_vs_numerical_gradient(testCase, 'U3',...
    %    @tensorsldaobj_matrixdata)
end


function test_parafac_lda_gradient_nr(testCase)    
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));

    analytical_vs_numerical_gradient(testCase, 'U1',...
        @parafacldaobj_matrixdata_normsratio)
    analytical_vs_numerical_gradient(testCase, 'U2',...
        @parafacldaobj_matrixdata_normsratio)
end

function analytical_vs_numerical_gradient(testCase, matrix_name, ...
    objective_function)
    p = 5;
    q = 7;
    r = 4;
    data_dimensions = [p, q];
    lower_dimensions = [3, 3];
    nsamples = 100;
    [Xs, ys] = get_data(nsamples, data_dimensions);
    
    [classmeandiffs, observationdiffs, nis] = ...
        classbased_differences(Xs, ys);
    
    myfun = @(U) objective_function(U,...
        classmeandiffs, observationdiffs, nis, lower_dimensions);

    initial_matrices_struct = get_initial_matrices(data_dimensions, lower_dimensions);
    numerically_estimated_gradient = get_numerically_estimated_gradient(...
        myfun, initial_matrices_struct, matrix_name);
    
    [F, analyticalderiv] = myfun(initial_matrices_struct);
    
    actual_difference = norm(...
        numerically_estimated_gradient-...
        analyticalderiv.(matrix_name), 'fro')/...
        norm(analyticalderiv.(matrix_name), 'fro');
    
    verifyEqual(testCase, actual_difference, 0, 'AbsTol', 1e-5)
end

function matrices_in_struct = get_initial_matrices(data_dimensions, Ks)
    for imode = 1:length(Ks)
        k = Ks(imode);
        p = data_dimensions(imode);

    [U, ~, ~] = svd(randn(p));
    matrices_in_struct.(['U', num2str(imode)]) = U(:, 1:k);


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

