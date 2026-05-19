#!/usr/bin/env python3
"""
Validation du profil de Blasius - Comparaison 3 methodes numeriques
Equation : f''' + 0.5*f*f'' = 0
CL : f(0)=0, f'(0)=0, f'(inf)=1

Methodes comparees :
  1. Runge-Kutta 4 (RK4)   - reference
  2. Euler explicite
  3. Runge-Kutta 2 (Heun)

Lecture des donnees depuis blasius_profile.dat genere par Oracle3D/usercod.f90

Usage : python3 blasius_validation.py
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# 1) RESOLUTION NUMERIQUE DE L EQUATION DE BLASIUS
# ============================================================

def rk4_step(f0, f1, f2, deta):
    """Un pas RK4 pour f'''=-0.5*f*f''"""
    k1f0 = deta * f1
    k1f1 = deta * f2
    k1f2 = deta * (-0.5 * f0 * f2)

    k2f0 = deta * (f1 + 0.5*k1f1)
    k2f1 = deta * (f2 + 0.5*k1f2)
    k2f2 = deta * (-0.5*(f0+0.5*k1f0)*(f2+0.5*k1f2))

    k3f0 = deta * (f1 + 0.5*k2f1)
    k3f1 = deta * (f2 + 0.5*k2f2)
    k3f2 = deta * (-0.5*(f0+0.5*k2f0)*(f2+0.5*k2f2))

    k4f0 = deta * (f1 + k3f1)
    k4f1 = deta * (f2 + k3f2)
    k4f2 = deta * (-0.5*(f0+k3f0)*(f2+k3f2))

    f0_new = f0 + (k1f0 + 2*k2f0 + 2*k3f0 + k4f0) / 6.0
    f1_new = f1 + (k1f1 + 2*k2f1 + 2*k3f1 + k4f1) / 6.0
    f2_new = f2 + (k1f2 + 2*k2f2 + 2*k3f2 + k4f2) / 6.0
    return f0_new, f1_new, f2_new


def euler_step(f0, f1, f2, deta):
    """Un pas Euler explicite"""
    f0_new = f0 + deta * f1
    f1_new = f1 + deta * f2
    f2_new = f2 + deta * (-0.5 * f0 * f2)
    return f0_new, f1_new, f2_new


def heun_step(f0, f1, f2, deta):
    """Un pas Runge-Kutta 2 (Heun)"""
    # Predictor (Euler)
    f0p = f0 + deta * f1
    f1p = f1 + deta * f2
    f2p = f2 + deta * (-0.5 * f0 * f2)
    # Corrector
    f0_new = f0 + 0.5*deta*(f1  + f1p)
    f1_new = f1 + 0.5*deta*(f2  + f2p)
    f2_new = f2 + 0.5*deta*((-0.5*f0*f2) + (-0.5*f0p*f2p))
    return f0_new, f1_new, f2_new


def shooting(step_func, f2_init, N, eta_max, tol=1e-10, max_iter=100):
    """Methode de tir pour trouver f''(0) tel que f'(inf)=1"""
    deta = eta_max / N
    f2_low, f2_high = 0.0, 2.0

    for it in range(max_iter):
        f2_init = 0.5*(f2_low + f2_high)
        f0, f1, f2 = 0.0, 0.0, f2_init
        for _ in range(N):
            f0, f1, f2 = step_func(f0, f1, f2, deta)
        if abs(f1 - 1.0) < tol:
            break
        if f1 < 1.0:
            f2_low  = f2_init
        else:
            f2_high = f2_init

    return f2_init, it+1


def integrate_blasius(step_func, N, eta_max):
    """Integre l equation de Blasius et retourne eta, f, f', f'' """
    f2_init, niter = shooting(step_func, 0.332, N, eta_max)
    deta = eta_max / N
    etas = np.zeros(N+1)
    F    = np.zeros(N+1)  # f
    FP   = np.zeros(N+1)  # f'  = u/U_inf
    FPP  = np.zeros(N+1)  # f''

    f0, f1, f2 = 0.0, 0.0, f2_init
    F[0], FP[0], FPP[0] = f0, f1, f2

    for i in range(N):
        f0, f1, f2 = step_func(f0, f1, f2, deta)
        etas[i+1] = (i+1)*deta
        F[i+1]    = f0
        FP[i+1]   = f1
        FPP[i+1]  = f2

    return etas, F, FP, FPP, f2_init, niter


# ============================================================
# 2) PARAMETRES
# ============================================================

N       = 2000       # nombre de pas d integration
ETA_MAX = 10.0       # eta_max (f' -> 1 pour eta > 5)

print("=" * 60)
print("  VALIDATION DU PROFIL DE BLASIUS")
print("  Equation : f''' + 0.5*f*f'' = 0")
print("=" * 60)

# ============================================================
# 3) INTEGRATION PAR LES 3 METHODES
# ============================================================

methods = [
    ("RK4   (Runge-Kutta 4)", rk4_step,   '-',  'blue'),
    ("RK2   (Heun)          ", heun_step,  '--', 'green'),
    ("Euler (ordre 1)       ", euler_step, ':',  'red'),
]

results = {}
for name, func, ls, col in methods:
    eta, F, FP, FPP, f2_0, nit = integrate_blasius(func, N, ETA_MAX)
    results[name] = (eta, F, FP, FPP, f2_0)
    print(f"  {name} : f''(0) = {f2_0:.8f}  ({nit} iterations shooting)")

print(f"\n  Valeur analytique  : f''(0) = 0.33205734 (Howarth 1938)")

# ============================================================
# 4) LECTURE DU FICHIER ORACLE3D (si present)
# ============================================================

oracle_data = None
if os.path.exists('blasius_profile.dat'):
    print("\n  Lecture de blasius_profile.dat (Oracle3D) ...")
    try:
        oracle_data = np.loadtxt('blasius_profile.dat', comments='#')
        print(f"  => {len(oracle_data)} points lus")
    except Exception as e:
        print(f"  => Erreur : {e}")

# ============================================================
# 5) SAUVEGARDE DES DONNEES
# ============================================================

# Sauvegarder les 3 profils dans un fichier
with open('blasius_comparison.dat', 'w') as f:
    f.write("# Profil de Blasius - Comparaison 3 methodes\n")
    f.write("# Equation : f''' + 0.5*f*f'' = 0\n")
    f.write("# CL : f(0)=f'(0)=0, f'(inf)=1\n")
    f.write("# Parametre de similarite : eta = psi * sqrt(U_inf / (nu * xi))\n")
    f.write("# Vitesse normalisee : u/U_inf = f'(eta)\n")
    f.write("#\n")
    f.write("# Colonnes : eta  f'_RK4  f'_RK2  f'_Euler  f_RK4  f''_RK4\n")
    f.write("#\n")
    eta_ref = results[list(results.keys())[0]][0]
    fp_rk4  = results[list(results.keys())[0]][2]
    fp_rk2  = results[list(results.keys())[1]][2]
    fp_eul  = results[list(results.keys())[2]][2]
    f_rk4   = results[list(results.keys())[0]][1]
    fpp_rk4 = results[list(results.keys())[0]][3]
    for i in range(len(eta_ref)):
        f.write(f"{eta_ref[i]:.6f}  {fp_rk4[i]:.8f}  {fp_rk2[i]:.8f}  {fp_eul[i]:.8f}  "
                f"{f_rk4[i]:.8f}  {fpp_rk4[i]:.8f}\n")

print(f"\n  Donnees sauvegardees dans : blasius_comparison.dat")
print(f"  Colonnes : eta | f'_RK4 | f'_RK2 | f'_Euler | f_RK4 | f''_RK4")

# ============================================================
# 6) TRACÉS
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 6))
fig.suptitle("Profil de Blasius — Comparaison des méthodes numériques\n"
             r"$f''' + \frac{1}{2}ff'' = 0$,  $f(0)=f'(0)=0$,  $f'(\infty)=1$",
             fontsize=13)

# --- Graphe 1 : u/U_inf = f'(eta) ---
ax = axes[0]
for name, func, ls, col in methods:
    eta, F, FP, FPP, f2_0 = results[name]
    label = f"{name.strip()} [f''(0)={f2_0:.5f}]"
    ax.plot(FP, eta, linestyle=ls, color=col, label=label, linewidth=1.8)

if oracle_data is not None:
    ax.plot(oracle_data[:,1], oracle_data[:,0], 'ko', markersize=3,
            label='Oracle3D', zorder=5)

ax.axvline(x=1.0, color='gray', linestyle=':', alpha=0.5)
ax.axhline(y=5.0, color='gray', linestyle=':', alpha=0.5, label=r'$\eta=5$ (99% U$_\infty$)')
ax.set_xlabel(r"$u/U_\infty = f'(\eta)$", fontsize=12)
ax.set_ylabel(r"$\eta = \psi\sqrt{U_\infty/(\nu\xi)}$", fontsize=12)
ax.set_title("Profil de vitesse longitudinale", fontsize=11)
ax.set_xlim(0, 1.1)
ax.set_ylim(0, ETA_MAX)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# --- Graphe 2 : Erreur relative par rapport a RK4 ---
ax = axes[1]
eta_ref = results[list(results.keys())[0]][0]
fp_ref  = results[list(results.keys())[0]][2]

for name, func, ls, col in methods[1:]:  # Euler et RK2 vs RK4
    eta, F, FP, FPP, f2_0 = results[name]
    # Erreur relative (eviter division par zero)
    err = np.abs(FP - fp_ref) / (np.abs(fp_ref) + 1e-10)
    ax.semilogy(err, eta, linestyle=ls, color=col,
                label=f"{name.strip()} vs RK4", linewidth=1.8)

ax.set_xlabel(r"Erreur relative $|f'_{method} - f'_{RK4}| / |f'_{RK4}|$", fontsize=11)
ax.set_ylabel(r"$\eta$", fontsize=12)
ax.set_title("Erreur relative / RK4 (référence)", fontsize=11)
ax.set_ylim(0, ETA_MAX)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# --- Graphe 3 : f, f', f'' avec RK4 ---
ax = axes[2]
eta, F, FP, FPP, f2_0 = results[list(results.keys())[0]]
ax.plot(F,   eta, 'b-',  label=r"$f(\eta)$",   linewidth=1.8)
ax.plot(FP,  eta, 'r-',  label=r"$f'(\eta) = u/U_\infty$", linewidth=1.8)
ax.plot(FPP, eta, 'g-',  label=r"$f''(\eta)$", linewidth=1.8)
ax.axhline(y=5.0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel(r"$f,\ f',\ f''$", fontsize=12)
ax.set_ylabel(r"$\eta$", fontsize=12)
ax.set_title("Fonctions de Blasius (RK4)", fontsize=11)
ax.set_ylim(0, ETA_MAX)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('blasius_validation.png', dpi=150, bbox_inches='tight')
print(f"  Figure sauvegardee : blasius_validation.png")
plt.show()

print("\n" + "="*60)
print("  RESUME")
print("="*60)
print(f"  f''(0) RK4   = {results[list(results.keys())[0]][4]:.8f}")
print(f"  f''(0) RK2   = {results[list(results.keys())[1]][4]:.8f}")
print(f"  f''(0) Euler = {results[list(results.keys())[2]][4]:.8f}")
print(f"  f''(0) ref   = 0.33205734  (Howarth 1938)")
print(f"  eta_99 (u=0.99*U_inf) ~ 5.0")