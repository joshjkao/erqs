import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from dotenv import load_dotenv


def run_command(cmd):
    return subprocess.run(cmd, check=False)


if __name__ == "__main__":
    load_dotenv()

    h_file = "instances/cluster_chain_n4_e-4.0000.txt"

    scratch_dir = os.getenv("SCRATCH_DIR", ".")

    # lambdas = np.linspace(0.00, 1.0, 8)
    lambdas = [1.0]
    n_updates = 100
    start_temp = 10.0
    end_temp = 10.0
    n_trials = 10

    data_dir = Path(scratch_dir) / Path(h_file).stem
    if not Path(data_dir).exists():
        Path(data_dir).mkdir(parents=True)

    commands = []
    instances = []

    for l in lambdas:
        for trial in range(n_trials):
            out_file = data_dir / f"out_l{l}_t{trial}.csv"
            instance = {
                "lambda": l,
                "n_updates": n_updates,
                "h_file": h_file,
                "start_temp": start_temp,
                "end_temp": end_temp,
                "trial": trial,
                "out_file": out_file,
            }
            instances.append(instance)

            commands.append(
                [
                    "../build/opt",
                    "--hamiltonian",
                    h_file,
                    "--updates",
                    f"{n_updates}",
                    "--lambda",
                    f"{l}",
                    "--output",
                    f"{out_file}",
                    "--start-temperature",
                    f"{start_temp}",
                    "--end-temperature",
                    f"{end_temp}",
                ]
            )

    with open("sweep_tasks.txt", "w") as f:
        f.writelines([" ".join(c) for c in commands])

    instances_df = pd.DataFrame(instances)
    instances_df.to_csv(Path(data_dir) / "instances.csv")

    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(run_command, commands))
