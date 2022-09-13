import pandas as pd

df = pd.read_csv("eso_lcs.csv")
print(df.set_index(["target_name", "fileid_img"]))
