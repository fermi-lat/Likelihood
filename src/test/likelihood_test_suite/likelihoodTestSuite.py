import os
import ctypes
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
import tracemalloc

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

libc = ctypes.CDLL(None)

# Line buffer Python stdout even when not writing to a terminal so output
# interleaves the same way it does on a terminal.
sys.stdout.reconfigure(line_buffering=True)

def flush_all():
    # Flush Python and C stdio buffers so output lands in the right place when the
    # file descriptors are switched.  When stdout is not a terminal these are block
    # buffered.  The C flush also covers Likelihood's C++ output and Minuit (f2c).
    sys.stdout.flush()
    sys.stderr.flush()
    libc.fflush(None)

# set up some directories for the functions to use
dirs = {}
dirs["FERMI_DIR"] = os.environ['FERMI_DIR']
dirs["GLOBAL_DATA_DIR"] = dirs['FERMI_DIR'] + "/refdata/fermi/"
dirs["LOCAL_DATA"] = os.environ['CONDA_PREFIX'] + "/share/fermitools/data/Likelihood/outref/"
key_path = "key_data"

print("\nPerforming time calibration ... ", end="")
scale = time_calib.benchmark()/time_calib.NOMINAL
# Experimental: more detailed calibration output, left disabled.
# print(f"scaling factor: {time_calib.benchmark():.3f}/{time_calib.NOMINAL:.3f} = {scale:.3f}")
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
        flush_all()
        os.dup2(outfile.fileno(),sys.stdout.fileno())
        os.dup2(outfile.fileno(),sys.stderr.fileno())
        # with redirect_stdout(f), redirect_stderr(f):
        start_time = time.time()
        try:
            # Experimental: tracemalloc (and the commented-out maxMemory line below) is
            # being evaluated as an alternative memory measure.  Its results are not
            # currently used.
            tracemalloc.start()
            tracemalloc.reset_peak()
            mem_usage = memory_usage(lambda:eval(cmd), backend="psutil_uss")
            current,peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            # print(mem_usage)
        except BaseException as e:
            print(''.join(traceback.format_tb(e.__traceback__)))
            print(e)
            test_ran = False
        stop_time = time.time()
        flush_all()
        os.dup2(stdout_fd,sys.stdout.fileno())
        os.dup2(stderr_fd,sys.stderr.fileno())

    # The data location varies by install, so replace it with the relative "data/" form
    # used in the key files. This keeps keys created from .out files portable.
    with open(f"output/{test}.out",'r') as file:
        output = file.read()
    with open(f"output/{test}.out",'w') as file:
        file.write(output.replace(dirs["LOCAL_DATA"], "data/"))

    if test_ran:
        if os.path.exists(f"{key_path}/{test}.key"):
            with open(f"output/{test}.out",'r') as file:
                output = file.read()
            with open(f"{key_path}/{test}.key",'r') as file:
                key = file.read()
            diff_results = "PASS" if key == output else "FAIL"

            elapsed = stop_time-start_time
            nominalTime = eval(f"mod.{test.split('.')[0]}.NOMINAL_TIME") * scale
            nominalMemory = eval(f"mod.{test.split('.')[0]}.NOMINAL_MEMORY")
            maxMemory = max(mem_usage[1:]) # The first point always seems high. The second value is always a drop to near zero and then it climbs again.
            # maxMemory = peak/1024/6 
            results[test]={
                'output_check': diff_results,
                'elapsed_time': elapsed,
                'max_memory': maxMemory,
                'nominal_time': nominalTime,
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
