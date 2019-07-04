import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.TestRunner
import matlab.unittest.TestSuite;

suite = matlab.unittest.TestSuite.fromName(...
'test_classification/test_cmda_discrimination');


runner = TestRunner.withTextOutput;
%runner.addPlugin(StopOnFailuresPlugin)
result = runner.run(suite)
%result = runner.run(test_classification, 'test_parafac_lda_discrimination')