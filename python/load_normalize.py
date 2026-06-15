import pandas as pd
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


if __name__ == "__main__":
    df = load_normalize("out")
    print(df.head())
