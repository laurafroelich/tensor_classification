function tests = test_get_diagonal_indices
    tests = functiontests(localfunctions);
end


function test_diagonal_indices_2d_dim3(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = magic(3);
    diag_inds = get_diagonal_indices(3, 2);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [8, 5, 2]'))
end


function test_diagonal_indices_2d_dim5(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = magic(5);
    diag_inds = get_diagonal_indices(5, 2);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [17, 5, 13, 21, 9]'))
end