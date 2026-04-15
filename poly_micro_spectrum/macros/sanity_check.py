import csv

f_sum = 0.0
d_sum = 0.0

with open("../output/poly_spectrum.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        f_sum += float(row["f_y"])
        d_sum += float(row["d_y"])

print("sum(f_y) =", f_sum)
print("sum(d_y) =", d_sum)