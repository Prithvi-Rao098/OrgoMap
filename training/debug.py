import csv
with open(CSV_PATH = "mechanism.csv", newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    print("CSV fieldnames:", reader.fieldnames)  # <--- Add this line to debug


