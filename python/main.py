import subprocess
import json
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from itertools import product


def run_benchark(max_depth, n_terms, n_factors, n_bits):
    # Construct the command line arguments
    # Ensure all arguments are strings
    cmd = [
        "../benchmark",
        str(max_depth),
        str(n_terms),
        str(n_factors),
        str(n_bits),
    ]

    # Run the process and capture stdout
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check for errors in the C++ execution
    if result.returncode != 0:
        return {"error": f"Process failed: {result.stderr}", "params": cmd}

    # Parse the JSON output from stdout
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        return {"error": "Invalid JSON output", "raw": result.stdout}


# 1. Helper function to handle a single task
def run_single_benchmark(args):
    # Unpack the parameters
    max_depth, n_factors, n_terms, n_bits = args
    # Call your existing C++ wrapper
    return run_benchark(max_depth, n_terms, n_factors, n_bits)


def run_all_parallel():
    n_bits = 10
    n_trials = 1000

    # 2. Generate all combinations of parameters (The "Flattening")
    # This creates a generator of tuples: (max_depth, n_factors, n_terms, n_bits)
    param_combinations = product(
        range(3, 8),  # max_depth
        range(3, 8),  # n_factors
        range(3, 8),  # n_terms
        [n_bits],
    )

    # Repeat the combinations for the number of trials
    tasks = [params for params in param_combinations for _ in range(n_trials)]

    # 3. Use ProcessPoolExecutor to run tasks in parallel
    # max_workers defaults to the number of processors on your machine
    print(f"Starting {len(tasks)} benchmarks across all available cores...")

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(run_single_benchmark, tasks))

    # 4. Save and return
    df = pd.DataFrame(results)
    df.to_csv("out.csv", index=False)
    return df


def run_all():
    results = []
    n_bits = 10
    n_trials = 1000
    for max_depth in range(3, 7):
        for n_factors in range(3, 7):
            for n_terms in range(3, 7):
                for trial in range(n_trials):
                    print(
                        f"trial {trial} for max_depth: {max_depth}, n_terms: {n_terms}, n_factors: {n_factors}"
                    )
                    result = run_benchark(max_depth, n_terms, n_factors, n_bits)
                    results.append(result)

    df = pd.DataFrame(results)
    df.to_csv("out")
    return df


def main():
    # run_all()
    run_all_parallel()


if __name__ == "__main__":
    main()
