function tests = test_classification_nd
    tests = functiontests(localfunctions);
end

function [Xs, ys] = get_simple_data(data_dimensions)
    nsamples = 200;
    [Xs, ys] = get_data(nsamples, data_dimensions);
end

function project_matrices_verify_predictions(testCase, ...
    projection_learning_function, prediction_function, data_dimensions, ...
    lower_dimensions)

    [Xs, ys] = get_simple_data(data_dimensions);
    Xs_test = Xs;
    ys_test = ys;
    
    Us = projection_learning_function(Xs, ys);
    
    predicted_probabilities = prediction_function(Xs, Xs_test, ys,...
        ys_test, lower_dimensions, Us);
    
    [~, predictions] = max(predicted_probabilities, [], 2);
    
    accuracy = mean(predictions == ys);
    
    verifyEqual(testCase, accuracy, 1, 'AbsTol', 1e-1)
end


%%%%%%%%%%%%%%%%%%
% Tests for the projection methods
%%%%%%%%%%%%%%%%%%

function test_cmda_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [4, 1, 2, 3];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) CMDA(Xs, ys, [], lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_parafac_lda_discrimination(testCase)
    import matlab.unittest.plugins.StopOnFailuresPlugin
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/manopt/', 'IncludeSubfolders', true));

    parafac_structure = true;
    
    data_dimensions = [5, 7, 3];
    lower_dimensions = [2, 2, 2];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) ManPDA(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_parafac_norms_ratio_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/manopt/', 'IncludeSubfolders', true));

    parafac_structure = true;
    
    data_dimensions = [5, 7];
    lower_dimensions = [2, 2];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) ManPDA_normsratio(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_tucker_lda_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/manopt/', 'IncludeSubfolders', true));

    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [4, 1, 2, 3];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) ManTDA(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_tucker_norms_ratio_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/manopt/', 'IncludeSubfolders', true));

    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [4, 1, 2, 3];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) ManTDA_normsratio(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_dgtda_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [2, 2, 2, 2];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) DGTDA(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_dater_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [4, 1, 2, 3];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) DATER(Xs, ys, [], lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_datereig_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [4, 1, 2, 3];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) DATEReig(Xs, ys, [], lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end


function test_hoda_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    
    data_dimensions = [5, 7, 3, 5];
    lower_dimensions = [2, 2, 2, 2];
    
    project_matrices_verify_predictions(testCase, ...
    @(Xs, ys) HODA(Xs, ys, lower_dimensions), ...
    @(Xs, Xs_test, ys, ys_test, k, Us) ...
    project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure),...
    data_dimensions, lower_dimensions)
end

