function tests = test_classification_nd
    tests = functiontests(localfunctions);
end

function test_cmda_discrimination(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    testCase.applyFixture(PathFixture('../../../../matlab_additions/02582nway_models/', 'IncludeSubfolders', true));
    
    parafac_structure = false;
    k = 2;
    
    [Xs, ys] = get_simple_data();
    
    CMDA(Xs, ys, [], [k, k, k])

    
    %project_matrices_verify_predictions(testCase, ...
    %@(Xs, ys) CMDA(Xs, ys, [], [k, k]), ...
    %@(Xs, Xs_test, ys, ys_test, k, Us) ...
    %project_and_predict(Xs, Xs_test, ys, ys_test, k, Us, parafac_structure))
end


function [Xs, ys] = get_simple_data()
    p = 5;
    q = 7;
    r = 3;
    nsamples = 200;
    [Xs, ys] = get_data(nsamples, [p, q, r]);
end

function project_matrices_verify_predictions(testCase, ...
    projection_learning_function, prediction_function, varargin)

    [Xs, ys] = get_simple_data();
    Xs_test = Xs;
    ys_test = ys;
    k = 2;
    
    Us = projection_learning_function(Xs, ys);
    
    predicted_probabilities = prediction_function(Xs, Xs_test, ys,...
        ys_test, k, Us);
    
    [~, predictions] = max(predicted_probabilities, [], 2);
    
    accuracy = mean(predictions == ys);
    
    verifyEqual(testCase, accuracy, 1, 'AbsTol', 1e-1)
end