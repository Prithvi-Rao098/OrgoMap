from kaggle.api.kaggle_api_extended import KaggleApi
import pandas as pd
import time

api = KaggleApi()
api.authenticate()

print("Starting download...")
start_time = time.time()

try:
    api.dataset_download_files(
        "broccoli201/chemical-reactions",
        path="./data",
        unzip=True,
        quiet=False  # Show progress
    )
    print(f"Download completed in {time.time()-start_time:.2f} seconds")
    
    df = pd.read_csv("./data/chemical-reactions.csv")
    print("Data loaded successfully:")
    print(df.head())
    
except Exception as e:
    print(f"Failed after {time.time()-start_time:.2f} seconds")
    print("Error:", str(e))