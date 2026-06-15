import pandas as pd
import numpy as np
import json


def load_normalize(path):
    # Assuming your json string is stored in a variable called json_string
    with open(path, "r") as f:
        data = json.load(f)

    # For this example, let's load it from a string
    # data = json.loads(json_string)

    # 1. Strip the outer list
    trials = data[0]

    # 2. Flatten the nested structure
    df = pd.json_normalize(
        trials,
        record_path=["fast_metrics"],
        meta=[
            "time_slow",
            "imag_inner",
            "real_inner",
            ["args", "h_terms"],
            ["args", "max_depth"],
            ["args", "n_factors"],
            ["args", "n_terms"],
        ],
        record_prefix="fast_",  # This prevents column name collisions!
    )

    return df


def quick_metrics(path):
    df = load_normalize(path)

    print(df.head())

    print(f"max real error: {np.max(df['fast_real_error'])}")
    print(f"max imag error: {np.max(df['fast_imag_error'])}")

    df_naive = df[df["fast_policy"] == "naive"]
    df_random = df[df["fast_policy"] == "random"]
    df_overlap = df[df["fast_policy"] == "overlap"]

    print(f"mean time for naive order: {np.mean(df_naive['fast_time_fast'])}")
    print(f"mean time for random order: {np.mean(df_random['fast_time_fast'])}")
    print(f"mean time for overlap order: {np.mean(df_overlap['fast_time_fast'])}")


if __name__ == "__main__":
    df = load_normalize("out")
    print(df.head())
