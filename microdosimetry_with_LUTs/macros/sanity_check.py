import csv, math

y = []
fy = []
yd = []

with open("../output/poly_spectrum.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        y.append(float(row["y_keV_per_um"]))
        fy.append(float(row["f_y"]))
        yd.append(float(row["yd_y"]))

# reconstruct C from the same 100-edge grid used in the code
edges = [10**(math.log10(0.1) + i*(math.log10(1000)-math.log10(0.1))/99) for i in range(100)]
widths = [edges[i+1] - edges[i] for i in range(99)]
C = math.log(10.0) * (math.log10(widths[1]) - math.log10(widths[0]))

print("C * sum(y * f(y))   =", C * sum(yi * fi for yi, fi in zip(y, fy)))
print("C * sum(y * d(y))   =", C * sum(ydi for ydi in yd))