function tests = test_cell_array_to_nd
    tests = functiontests(localfunctions);
end


function test_cell_array_to_nd_array_arbitraryn_modes(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    n_obs = 5;
    n_modes = randsample(2:10, 1, true);
    obs_dimensions = randsample(1:8, n_modes, true);
    xs_nd_array = zeros([obs_dimensions, n_obs]);
    xs_cell = cell(n_obs);
    S.type = '()';
    subs = repmat({':'}, 1, n_modes);
    
    for iobs = 1:n_obs
        S.subs = [subs, iobs];
        rand_obs = normrnd(0, 1, obs_dimensions);
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
