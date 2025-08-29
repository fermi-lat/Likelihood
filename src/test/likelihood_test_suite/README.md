# Likelihood Test Suite
This tool is a regression/perfromance test suite for the likelihood tool.  The initial tests are all drawn from test scripts that Jean Ballet sent in while we were working on the memory use issues in the likelihood tool.  It is designed to be run against the current version of the tools during development to verify that changes have not impacted any of the know uses of the software.  It is extensible and additional tests can be added to the suite as they are developed.

The suite verifies three items:
* the correct output values are computed
* the execution time has not changed
* the memory used has not changed

Each of these are reported for each test run.

## Requirements
The script is written in Python and requres Python 3 to be installed. The only non-standard Python library required is the `memory_profiler` package that can be installed via pip or conda.

## Contents
The top level directory contains three files and four directories:
* **README.md** - this file
* **likelihoodTestSuite.py** - the main program file for the test suite 
* **time_calib.py** - a module containing a small function used for time calibration
* **data/** - directory containing the data used by the test scripts
* **key_data/** - directory containing the expected output files (keys) to verify correct execution of the various tests.
* **output/** - directory containing the output results for each script.  Contents are compared against the corresponding files in the key_data directory to validate test execution.
* **scripts/** - directory containing the individual tests to be run
  * **scripts/inactive/** - scripts in this directory are not executed.  If you wish to temporarily disable execution of a script, move it here.  It also contains a template file for writing new scripts
  
## Running the Test Suite
The suite is invoked by running:

```sh
python likelihoodTestSuite.py
```
No other input is needed.  

The program starts by running a short time calibration so that the local execution times can be compared to the nominal time in the scripts.

After that the program reads all the scripts in the `scripts/` directory and executes them, printing out the results of each test.  The output of the program looks like the following:

```text
Performing time calibration ... scaling factor: 0.980

|-------------------------------------|------|-------|-------|------------------|----------|----------|--------------------|
|                                     |Output|Elapsed|Nominal|                  |   Max    | Nominal  |                    |
| Test                                |Check |Time(s)|Time(s)|    Time Delta    |Memory(MB)|Memory(MB)|    Memory Delta    |
|-------------------------------------|------|-------|-------|------------------|----------|----------|--------------------|
| test_EBL                            | PASS |  3.40 |  3.14 |   0.20   (6.43%) |   116.02 |   115.00 |     1.02   (0.88%) |
| test_FT_border                      | PASS |  6.75 |  6.47 |   0.15   (2.38%) |  1114.76 |  1150.00 |   -35.24  (-3.06%) |
| test_FT_ub                          | PASS |  6.92 |  6.37 |   0.42   (6.55%) |  1216.09 |  1200.00 |    16.09   (1.34%) |
| test_UpperLimit                     | PASS |  0.49 |  0.49 |  -0.01  (-2.07%) |    86.16 |    87.00 |    -0.84  (-0.96%) |
| test_delete_src                     | PASS |  0.45 |  0.49 |  -0.05 (-10.17%) |    85.46 |    86.00 |    -0.54  (-0.63%) |
| test_eric                           | PASS |  0.56 |  0.64 |  -0.09 (-13.69%) |    87.65 |    88.00 |    -0.35  (-0.40%) |
| test_freeze_thaw                    | PASS |  2.92 |  2.94 |  -0.08  (-2.76%) |   105.07 |   106.00 |    -0.93  (-0.88%) |
| test_freeze_thaw_weighted           | PASS |  7.54 |  7.16 |   0.24   (3.42%) |   359.69 |   360.00 |    -0.31  (-0.09%) |
| test_freeze_thaw_weighted_clearmaps | PASS |  2.74 |  2.60 |   0.09   (3.29%) |   107.57 |   108.00 |    -0.43  (-0.40%) |
| test_memory                         | PASS | 19.21 | 20.59 |  -1.79  (-8.70%) |   895.39 |   895.00 |     0.39   (0.04%) |
| test_memory_tHb1                    | PASS |  0.80 |  0.83 |  -0.05  (-5.63%) |    94.73 |    95.00 |    -0.27  (-0.28%) |
| test_memory_weighted                | PASS | 25.47 | 27.45 |  -2.53  (-9.23%) |  1264.43 |  1265.00 |    -0.57  (-0.05%) |
| test_memory_weighted_clearmaps      | PASS |105.02 |105.87 |  -2.98  (-2.82%) |  1109.53 |  1110.00 |    -0.47  (-0.04%) |
| test_restore                        | PASS |  1.18 |  1.23 |  -0.07  (-5.62%) |   117.73 |   120.00 |    -2.27  (-1.89%) |
|-------------------------------------|------|-------|-------|------------------|----------|----------|--------------------|
```
The first line gives the results of the timing calibration.  This is the time scaling factor for the computer the program is being run on.  All nominal times in the test scripts are multiplied by this factor to get the nominal times expected for this run of the program.

The output columns in the results table are as follows:
* Test - the name of the executed tests
* Ouput Check - results of comparing the output to the key_data file for the tests.  PASS = correct output, FAIL = bad output
* Elapsed Time(s) - the time in seconds the test actually took to run
* Nominal Time(s) - the expected nominal time (in seconds) for the test, scaled by the computed scaling factor
* Time Delta - The difference between the elapsed and nominal time.  The difference is given in seconds as well as the percentage variation (in parentheses).
* Max Memory(MB) - the maximum memory used in execution of the tests, in Megabytes
* Nominal Memory(MB) - the expected memory usage in Megabytes
* Memory Delta - the difference between the actual and expected memory use. The difference is given in MB as well as the percentage variation (in parentheses).

The actual output will be color coded:
* For the `Output Check` column:
  * Green - PASS
  * Red - FAIL
* For the `Time Delta` and `Memory Delta` colums:
  * Green - actual execution value is within 10% of the nominal value
  * Red - actual value is more than 10% higher than the nominal value
  * Amber/Yellow - actual value is more than 10% lower than the nominal value

***Note:** The percentage for triggering the color change can be changed by setting the threshold parameter in the printing function calls. 10% is the default.*

## Adding a Test
To add a new test, you need to create a new test script file, create a key_file, and set the nominal time and memory values for the test.

***Note:** While creating a new test, you might want to move all the existing tests to the `scripts/inactive/` directory so you don't have to wait for them to all run while testing your new script.*

### Create the test script file
The main program is expecting to find files with names of the form `test_<testfilename>.py` in the `scripts/` directory.  Inside of those files it looks for a function with a name of the form `test_<testname>`.  The function should take one parameter which is a dictionary that contains directory information so the test can find the data directory.

There is a template test file in the `scripts/inactive/` directory that you can copy and modify.  The contents of that file are listed below:

```python
NOMINAL_MEMORY = 1
NOMINAL_TIME = 1

def test_template(dirs):
    # test data is in dirs['LOCAL_DATA']

    #Test code goes here
    pass
```
The `NOMINAL_MEMORY` and `NOMINAL_TIME` constants are where you set these values for this test.  They are given in MB and seconds respectively.

***Note:** The main program actually looks for any number of tests within a single test script file.  However, as currently written, it assumes the same nominal execution time and memory usage for all tests in the file so you should only put one test in each test file at the moment.*

Add whatever test code you want to execute for the test in the function (renaming it from `test_template` to something that describes the test).  Also be sure to rename the file from `test_template.py` to something reasonable.  The test filenames and test names should be unique.

### Create the Key file
The key file is the expected output from runnin the test.  The test suite caputres `STDOUT` and `STDERR` to a single file for comparison.  The key files are stored in the `key_data/` directory and have filenames of the form: `test_<testfilename>.test_<testname>.key` where `test_<testfilename>` is the name of the test script file and `test_<testname>` is the name of the test function in that file.

The easiest way to create the key file is to run the test suite in an environment with a version of the Fermitools that produces the correct results.  The program will print the message `**key file missing` for the test.  However, it will have generated the correct output in the `output/` directory.  In that directory, the files are named `test_<testfilename>.test_<testname>.out` simply change the `.out` extention to `.key` and move the file to the `key_data/` directory and you are good to go.

### Setting Nominal Time and Memory values
Now that you have a proper key file, run the test suite again.  This time it will print out the values for the `Elapsed Time` and `Max Memory` columns for your new test.  Set the `NOMINAL_MEMORY` and `NOMINAL_TIME` constants in your test script to the measured values.  You might want to run the test a few times and use the average of a few runs.

**Important:** When setting the `NOMINAL_TIME` value, remember that this is scaled by the scaling factor computed at the start of the run.  You need to **divide** your measured `Elapsed Time` by this scale factor to set the `NOMINAL_TIME` for your script.

Now you should be able to run the test suite and see green values for you new test.  You are ready to run it against development builds to verify that all is still well after changes.  Don't forget to move all the tests out of the `inactive/` folder if you moved them there for test script development.

## TODO
Possible updates to make to the program.
* Add command-line option to set the color coding threshold
* Make the nominal time and memory correspond to individual tests instead of test files so multiple related tests can be put in a single test file.
* Add option for results to be written to a file.