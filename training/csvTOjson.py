import pandas as pd
import json

df = pd.read_csv("your_chemistry_data.csv")

with open("converted_dataset.jsonl", "w") as f:
    for _, row in df.iterrows():
        example = {
            "contents": [
                {"role": "user", "parts": [{"text": row["prompt"]}]},
                {"role": "model", "parts": [{"text": row["response"]}]}
            ]
        }
        f.write(json.dumps(example) + "\n")
