function tests = test_classification_nd
    tests = functiontests(localfunctions);
end



function test_cell_array_to_nd_array_arbitraryn_modes(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    n_obs = 5;
    n_modes = randsample(2:4, 1);
    obs_dimensions = randsample(1:3, n_modes);
    xs_nd_array = zeros([obs_dimensions, n_obs]);
    xs_cell = cell(n_obs);
    S.type = '()';
    S.subs = [repmat({':'}, 1, n_modes)]
    
    for iobs = 1:n_obs
        rand_obs = normrnd(0, 1, obs_dimensions);
        %xs_nd_array(:,iobs) = rand_obs;
        xs_nd_array = subsasgn(xs_nd_array, S, rand_obs);
        xs_cell{iobs} = rand_obs;
        
    end
    
    result = cell_array_to_nd_array(xs_cell);
    
    isequal(result, xs_nd_array)
end

function test_cell_array_to_nd_array_2modes(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    n_obs = 5;
    obs_dimensions = [7, 3];
    xs_nd_array = zeros([obs_dimensions, n_obs]);
    xs_cell = cell(n_obs);
    
    for iobs = 1:n_obs
        rand_obs = normrnd(0, 1, obs_dimensions);
        xs_nd_array(:,:,iobs) = rand_obs;
        xs_cell{iobs} = rand_obs;
        
    end
    
    result = cell_array_to_nd_array(xs_cell);
    
    isequal(result, xs_nd_array)
end

function test_cell_array_to_nd_array_3modes(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    n_obs = 5;
    obs_dimensions = [7, 3, 10];
    xs_nd_array = zeros([obs_dimensions, n_obs]);
    xs_cell = cell(n_obs);
    
    for iobs = 1:n_obs
        rand_obs = normrnd(0, 1, obs_dimensions);
        xs_nd_array(:,:,:,iobs) = rand_obs;
        xs_cell{iobs} = rand_obs;
        
    end
    
    result = cell_array_to_nd_array(xs_cell);
    
    isequal(result, xs_nd_array)
end


function test_cell_array_to_nd_array_4modes(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    n_obs = 6;
    obs_dimensions = [2, 3, 4, 5];
    xs_nd_array = zeros([obs_dimensions, n_obs]);
    xs_cell = cell(n_obs);
    
    for iobs = 1:n_obs
        rand_obs = normrnd(0, 1, obs_dimensions);
        xs_nd_array(:,:,:, :, iobs) = rand_obs;
        xs_cell{iobs} = rand_obs;
        
    end
    
    result = cell_array_to_nd_array(xs_cell);
    
    isequal(result, xs_nd_array)
end

function wip_test_cmda_discrimination(testCase)
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