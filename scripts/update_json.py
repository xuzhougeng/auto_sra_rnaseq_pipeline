import os
import sys
import json


dict_key = sys.argv[1]
meta_file = sys.argv[2]

# update the json file
if os.path.exists("file_dict.json"):
    with open("file_dict.json", "r") as f:
        file_dict = json.load(f)
# get the key-value pair
else:
    file_dict = {}
    
file_dict[dict_key] = meta_file

print("Writing the processed file to file_dict.json") 
with open("file_dict.json", "w") as f:
    json.dump(file_dict, f)
