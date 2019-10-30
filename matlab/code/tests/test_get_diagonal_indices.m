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

function test_diagonal_indices_3d_dim3(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:27, [3, 3, 3]);
    diag_inds = get_diagonal_indices(3, 3);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 14, 27]'))
end


function test_diagonal_indices_4d_dim3(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:81, [3, 3, 3, 3]);
    diag_inds = get_diagonal_indices(3, 4);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 41, 81]'))
end

function test_diagonal_indices_4d_dim2(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:(2^4), [2, 2, 2, 2]);
    diag_inds = get_diagonal_indices(2, 4);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 16]'))
end



function test_diagonal_indices_6d_dim3(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:(3^6), [3, 3, 3, 3, 3, 3]);
    diag_inds = get_diagonal_indices(3, 6);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 365, 729]'))
end


function test_diagonal_indices_3d_dim5(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:(5^3), [5, 5, 5]);
    diag_inds = get_diagonal_indices(5, 3);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 32, 63, 94, 125]'))
end


function test_diagonal_indices_4d_dim7(testCase)
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('../', 'IncludeSubfolders', true));

    mat = reshape(1:(7^4), [7, 7, 7, 7]);
    diag_inds = get_diagonal_indices(7, 4);
    diag_elems = mat(diag_inds);
    
    assert(isequal(diag_elems, [1, 401, 801, 1201, 1601, 2001, 2401]'))
end
