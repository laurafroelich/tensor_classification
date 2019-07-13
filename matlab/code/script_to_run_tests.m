import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.TestRunner
import matlab.unittest.TestSuite;

suite_cmda = matlab.unittest.TestSuite.fromName(...
'test_classification/test_cmda_discrimination');

suite_datereig = matlab.unittest.TestSuite.fromName(...
'test_classification/test_datereig_discrimination');

suite_hoda = matlab.unittest.TestSuite.fromName(...
'test_classification/test_hoda_discrimination');

suite_dater = matlab.unittest.TestSuite.fromName(...
'test_classification/test_dater_discrimination');

suite_dgtda = matlab.unittest.TestSuite.fromName(...
'test_classification/test_dgtda_discrimination');

suite_bilinear_tucker = matlab.unittest.TestSuite.fromName(...
'test_classification/test_bilinear_tucker_discrimination');

suite_bilinear_parafac = matlab.unittest.TestSuite.fromName(...
'test_classification/test_bilinear_parafac_discrimination');

suite_lda_parafac = matlab.unittest.TestSuite.fromName(...
'test_classification/test_parafac_lda_discrimination');

suite_lda_parafac_nr = matlab.unittest.TestSuite.fromName(...
'test_classification/test_parafac_norms_ratio_discrimination');

runner = TestRunner.withTextOutput;
%runner.addPlugin(StopOnFailuresPlugin)

%result = runner.run([suite_datereig, suite_cmda, suite_hoda, ...
%    suite_dater, suite_dgtda, suite_bilinear_tucker, ...
%suite_bilinear_parafac, suite_lda_parafac, suite_lda_parafac_nr])

result = runner.run(suite_lda_parafac_nr)
