from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    results = Path("output/results/20260908103650")
    naccess_results = Path("output/naccess_results")

    with open(Path("data/pdb_ids.txt"), "r") as f:
        ids = [line.strip() for line in f.readlines()]

    resultats = {
        "custom" : {},
        "naccess" : {}
    }

    for identifiant in ids:
        with open(results / Path(f"{identifiant}/{identifiant}_rsa.txt"), "r") as f:
            resultats["custom"][identifiant] = float(f.readlines()[-1][6:])

        with open(naccess_results / Path(f"{identifiant}/.rsa"), "r") as f:
            resultats["naccess"][identifiant] = float(f.readlines()[-1][6:22])

    print(resultats)

    plt.figure(figsize=(7, 7))

    plt.scatter(list(resultats["custom"].values()), list(resultats["naccess"].values()))

    # droite y = x
    lims = [
        min(min(list(resultats["naccess"].values())), min(list(resultats["custom"].values()))),
        max(max(list(resultats["naccess"].values())), max(list(resultats["custom"].values())))
    ]
    plt.plot(lims, lims, linestyle="--")

    plt.xlabel("Surface NACCESS")
    plt.ylabel("Surface calculée")
    plt.title("Comparaison des surfaces")

    plt.xlim(lims)
    plt.ylim(lims)
    plt.axis("equal")
    plt.show()


    moyenne = (
        np.asarray(list(resultats["naccess"].values())) + np.asarray(list(resultats["custom"].values()))
    ) / 2

    difference = (
        np.asarray(list(resultats["custom"].values())) - np.asarray(list(resultats["naccess"].values()))
    )

    biais = difference.mean()
    sd = difference.std()

    plt.figure(figsize=(8, 5))

    plt.scatter(moyenne, difference)

    plt.axhline(biais, linestyle="--", label=f"Biais = {biais:.2f}")
    plt.axhline(biais + 1.96 * sd, linestyle=":", label="± 1.96 SD")
    plt.axhline(biais - 1.96 * sd, linestyle=":")

    plt.xlabel("Moyenne des deux méthodes")
    plt.ylabel("Différence (code − NACCESS)")
    plt.title("Bland–Altman : comparaison des surfaces")
    plt.legend()

    plt.show()
