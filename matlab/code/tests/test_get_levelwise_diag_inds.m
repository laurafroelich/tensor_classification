function tests = test_get_levelwise_diag_inds
    tests = functiontests(localfunctions);
end


function test_get_levelwise_diag_inds_2d_dim3_1level(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = magic(3);
    diag_inds =  get_level_wise_diag_inds(3, 1, 2);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [8, 5, 2]'))
end


function test_get_levelwise_diag_inds_2d_dim3_2levels(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));
    
    mat = zeros(3, 3, 2);
    mat(:,:,1) = magic(3);
    mat(:,:,2) = magic(3);
    diag_inds =  get_level_wise_diag_inds(3, 2, 2);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [8, 5, 2, 8, 5, 2]'))
end


function test_get_levelwise_diag_inds_3d_dim3_4levels(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = zeros(3, 3, 3, 4);
    for ilevel = 1:4
    mat(:, :, :, ilevel) = reshape(1:27, [3, 3, 3]);
    end
    diag_inds = get_level_wise_diag_inds(3, 4, 3);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 14, 27, 1, 14, 27, 1, 14, 27, 1, 14, 27]'))
end