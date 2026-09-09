from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


if __name__ == "__main__":
    molecule = "3I40"
    results = Path(f"output/results_test/{molecule}")
    naccess_results = Path(f"output/naccess_results/{molecule}")


    resultats = {
        "custom" : {},
        "naccess" : {}
    }

    with open(results / Path(f"{molecule}_rsa.txt"), "r") as f:
        for line in f:
            if line.startswith("RES"):
                resultats["custom"][int(line[9:14])] = float(line[14:23])

    with open(naccess_results / ".rsa", "r") as f:
        for line in f:
            if line.startswith("RES"):
                resultats["naccess"][int(line[9:14])] = float(line[14:23])

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

    plt.figure(figsize=(8, 5))
    plt.plot(resultats["custom"].keys(), resultats["custom"].values(), label = "Custom")
    plt.plot(resultats["naccess"].keys(), resultats["naccess"].values(), label = "Naccess")
    plt.xticks([1] + list(range(5, max(resultats["naccess"].keys()) + 1, 5)))
    plt.show()
