import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.TestRunner
import matlab.unittest.TestSuite;

suite_cmda = matlab.unittest.TestSuite.fromName(...
'test_classification/test_cmda_discrimination');
suite_datereig = matlab.unittest.TestSuite.fromName(...
'test_classification/test_datereig_discrimination');
suite_hoda = matlab.unittest.TestSuite.fromName(...
'test_classification/test_hoda_discrimination');

runner = TestRunner.withTextOutput;
%runner.addPlugin(StopOnFailuresPlugin)

result = runner.run([suite_datereig, suite_cmda, suite_hoda])