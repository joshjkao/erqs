# ERQS Public Repo

## Building 

From the root directory:
```
cmake -CMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build -j<N_CPU>
```

## Running

The `main` executable is a demo and can be changed to whatever. Use it to test new features, debug, etc.
```
./build/main
```

The `benchmark` executable constructs a random quantum state based on the command line parameters, simplifies it, normalizes it, and computes its inner product with itself.
```
./build/benchmark <MAX_DEPTH> <TERMS_PER_SUM> <FACTORS_PER_PRODUCT> <N_QUBITS>
```

Sample output:
```
{
"max_depth": 5,
"n_terms": 5,
"n_factors": 5,
"n_bits": 10,
"n_nodes": 308,
"c_time": 2666,
"mem_max": 26712,
"resulting_terms": 1,
"inner": 0.889827
}
```

## Python Directory

Construct the virtual environment with $uv$:
```
cd python
uv sync
```

Open the $main.py$ to see how the parameter sweep is defined. As is, `run_all()` will run a parameter sweep using one cpu. `run_all_parallel()` will run a parameter sweep using all available cpus.

Either function will output a file $out.csv$ containing the results of each trial in the parameter sweep.

Example usage (from the `python` directory)
```
uv run main.py
uv run python
import pandas as pd
import import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("out.csv")
sns.pointplot(data=df, x="n_nodes", y="c_time")
plt.show()
```
