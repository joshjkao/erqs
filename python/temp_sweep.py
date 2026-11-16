import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from dotenv import load_dotenv
from pathlib import Path


def run_command(cmd):
    return subprocess.run(cmd, check=False)


if __name__ == "__main__":
    load_dotenv()

    scratch_dir = os.getenv("SCRATCH_DIR", ".")

    lambdas = np.linspace(0.1, 2.0, 6)
    n_updates = 200
    h_file = "cluster_chain_n6_e-6.0000.txt"
    start_temp = 1.0
    end_temp = 0.1

    commands = []

    for l in lambdas:
        out_file = Path(scratch_dir) / f"out_l{l}.csv"
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
                f"{end_temp}"
            ]
        )

    with open("sweep_tasks.txt", "w") as f:
        f.writelines([c for c in commands])

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(run_command, commands))

    print("started thread pool")
