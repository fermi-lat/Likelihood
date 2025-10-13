import os
from contextlib import redirect_stdout, redirect_stderr
import re
import inspect
import time
import sys
import gc
import traceback
from memory_profiler import memory_usage
import time_calib
from pathlib import Path

#color codes
ANSI_RED = "\033[91m"
ANSI_GREEN = "\033[92m"
ANSI_AMBER = "\033[33m"
ANSI_CLEAR = "\033[0m"

def print_results_header():
    print()
    print_results_divider()
    print("|                                     |Output|Elapsed|Nominal|                  |   Max    | Nominal  |                    |")
    print("| Test                                |Check |Time(s)|Time(s)|    Time Delta    |Memory(MB)|Memory(MB)|    Memory Delta    |")
    print_results_divider()

def print_test_name(name):
    print(f"| {name.split('.')[0]:36s}|", end="", flush=True)

def change_output_color(percent, threshold = 10):
    if percent > threshold:
        print(ANSI_RED, end="")
    elif percent < -threshold:
        print(ANSI_AMBER, end="")
    else:
        print(ANSI_GREEN, end="")
    
def print_results_data(data, threshold = 10):
    # output check data
    print((ANSI_GREEN if data['output_check'] == "PASS" else ANSI_RED), end="")
    print(f" {data['output_check']:5s}" + ANSI_CLEAR, end="")

    # time data
    print(f"|{data['elapsed_time']:6.2f} |{data['nominal_time']:6.2f} |", end="")
    percent = data['time_delta']/data['nominal_time']*100
    change_output_color(percent, threshold)
    percentString = f"({percent:.2f}%)"
    print(f"{data['time_delta']:7.2f} {percentString:>9s} " + ANSI_CLEAR, end="")

    #memory data
    print(f"|{data['max_memory']:9.2f} |{data['nominal_memory']:9.2f} |", end="")
    percent = data['memory_delta']/data['nominal_memory']*100
    change_output_color(percent, threshold)
    percentString = f"({percent:.2f}%)"
    print(f" {data['memory_delta']:8.2f} {percentString:>9s}" + ANSI_CLEAR + " |")

def print_results_divider():
    print("|-------------------------------------|------|-------|-------|------------------|----------|----------|--------------------|")

def print_full_results(results, threshold = 10):
    print_results_header()
    for key in results:
        print_test_name(key)
        print_results_data(results[key], threshold)
    print_results_divider()


stdout_fd = os.dup(sys.stdout.fileno())
stderr_fd = os.dup(sys.stderr.fileno())

# set up some directories for the functions to use
dirs = {}
dirs["FERMI_DIR"] = os.environ['FERMI_DIR']
dirs["GLOBAL_DATA_DIR"] = dirs['FERMI_DIR'] + "/refdata/fermi/"
dirs["LOCAL_DATA"] = os.environ['CONDA_PREFIX'] + "/share/fermitools/data/test-scripts/"
key_path = os.environ['CONDA_PREFIX'] + "/share/fermitools/data/test-scripts/outref/"

print("\nPerforming time calibration ... ", end="")
scale = time_calib.benchmark()/time_calib.NOMINAL
print(f"scaling factor: {scale:.3f}")

# Create output directory if not present
OUTPUT_PATH = "output"
Path(OUTPUT_PATH).mkdir(parents=True, exist_ok=True)

# find all the tests in the scripts directory and build a list of the
# filenames without the .py extension
listing = os.listdir('scripts')
testFileRegEx = re.compile(r'(test.*)\.py')
testModules = []
for item in listing:
    match = testFileRegEx.search(item)
    if match:
        testModules.append(match.group(1))
# print(testModules)

# loop over all the test modules and import them.  They will all be
# sub elements of the mod object
for module in testModules:
    # print(module)
    mod = __import__(f"scripts.{module}")

# Now find all the actual tests in the loaded modules
testRegEx = re.compile(r'test.*')
tests = []
# Loop over all the testModules above
for name in testModules:
    element = getattr(mod, name)
    if inspect.ismodule(element):
        # print(f"Processing test module: {name}")
        contents = dir(element)
        for item in contents:
            if testRegEx.fullmatch(item):
                # print(f"Found test: {item}")
                tests.append(f"{name}.{item}")
# print(tests)

results = {}
print_results_header()
for test in sorted(tests):
    gc.collect()  # force garbage collection to make memory measurement more accurate
    # time.sleep(1)
    cmd = f"mod.{test}(dirs)"
    print_test_name(test)
    # print(f"\nRunning test: {test}:")
    test_ran = True
    with open(f"output/{test}.out", 'w') as outfile:
        os.dup2(outfile.fileno(),sys.stdout.fileno())
        os.dup2(outfile.fileno(),sys.stderr.fileno())
        # with redirect_stdout(f), redirect_stderr(f):
        start_time = time.time()
        try:
            mem_usage = memory_usage(lambda:eval(cmd), backend="psutil_uss")
            # print(mem_usage)
        except BaseException as e:
            print(''.join(traceback.format_tb(e.__traceback__)))
            print(e)
            test_ran = False
        stop_time = time.time()
        os.dup2(stdout_fd,sys.stdout.fileno())
        os.dup2(stderr_fd,sys.stderr.fileno())

    if test_ran:
        if os.path.exists(f"{key_path}/{test}.key"):
            with open(f"output/{test}.out",'r') as file:
                output = file.read()
            with open(f"{key_path}/{test}.key",'r') as file:
                key = file.read()
            diff_results = "PASS" if key == output else "FAIL"

            elapsed = stop_time-start_time
            nominalTime = eval(f"mod.{test.split('.')[0]}.NOMINAL_TIME")
            nominalMemory = eval(f"mod.{test.split('.')[0]}.NOMINAL_MEMORY")
            maxMemory = max(mem_usage[1:]) # The first point always seems high. The second value is always a drop to near zero and then it climbs again.
            results[test]={
                'output_check': diff_results,
                'elapsed_time': elapsed,
                'max_memory': maxMemory,
                'nominal_time': nominalTime * scale,
                'nominal_memory': nominalMemory,
                'time_delta': elapsed - nominalTime,
                'memory_delta': maxMemory - nominalMemory,
                'memory_profile': mem_usage
            }
            print_results_data(results[test])
        else:
            print(" **key file missing")
    else:
        print(" Test failed to run properly. See output for details.")

print_results_divider()
# print_full_results(results)
