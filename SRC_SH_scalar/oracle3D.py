import numpy as np
import matplotlib.pyplot as plt

filename = "field.plt"

# -------------------------------
# Lire en-tête
# -------------------------------

with open(filename, "r") as f:

    line = f.readline()

    # Lire VARIABLES
    while "VARIABLES" not in line.upper():
        line = f.readline()

    vars_line = line
    variables = vars_line.split('"')[1::2]

    print("Variables:", variables)

    # Lire ZONE
    line = f.readline()

    while "ZONE" not in line.upper():
        line = f.readline()

    zone_line = line

    print("Zone:", zone_line)

    # Extraire I, J, K
    import re

    I = int(re.search(r'I\s*=\s*(\d+)', zone_line).group(1))
    J = int(re.search(r'J\s*=\s*(\d+)', zone_line).group(1))
    K_match = re.search(r'K\s*=\s*(\d+)', zone_line)

    if K_match:
        K = int(K_match.group(1))
    else:
        K = 1

    print("Grid:", I, J, K)

    npts = I * J * K
    nvar = len(variables)

    # -------------------------------
    # Lire les données numériques
    # -------------------------------

    data = np.loadtxt(f)

# -------------------------------
# Reshape
# -------------------------------

data = data.reshape((npts, nvar))

# Dictionnaire des variables

fields = {}

for i, var in enumerate(variables):
    fields[var] = data[:, i].reshape((J, I))

# -------------------------------
# Extraire coordonnées
# -------------------------------

X = fields["X"]
Y = fields["Y"]

# -------------------------------
# POTENTIEL ELECTRIQUE VP
# -------------------------------

plt.figure()

plt.pcolormesh(X, Y, fields["VP"], shading="auto")

plt.colorbar(label="Potentiel VP")

plt.title("Potentiel électrique VP")

plt.xlabel("X")
plt.ylabel("Y")

plt.axis("equal")

plt.show()

# -------------------------------
# VITESSE U
# -------------------------------

plt.figure()

plt.pcolormesh(X, Y, fields["U"], shading="auto")

plt.colorbar(label="Vitesse U")

plt.title("Vitesse U")

plt.xlabel("X")
plt.ylabel("Y")

plt.axis("equal")

plt.show()

# -------------------------------
# NORME VITESSE
# -------------------------------

Umag = np.sqrt(
    fields["U"]**2 +
    fields["V"]**2 +
    fields["W"]**2
)

plt.figure()

plt.pcolormesh(X, Y, Umag, shading="auto")

plt.colorbar(label="|U|")

plt.title("Norme vitesse")

plt.xlabel("X")
plt.ylabel("Y")

plt.axis("equal")

plt.show()

# -------------------------------
# FORCE EHD
# -------------------------------

plt.figure()

plt.pcolormesh(X, Y, fields["Fmag"], shading="auto")

plt.colorbar(label="Force EHD")

plt.title("Force plasma")

plt.xlabel("X")
plt.ylabel("Y")

plt.axis("equal")

plt.show()