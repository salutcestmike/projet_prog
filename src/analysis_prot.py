from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr


if __name__ == "__main__":
    molecule = "168L"
    name = "Lysozyme"
    results = Path(f"output/results_test/{molecule}")
    naccess_results = Path(f"output/naccess_results/{molecule}")
    plt.rcParams.update({'font.size': 16})

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
    plt.title(f"Comparaison des surfaces, {name} {molecule}")

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
    print(sd)
    plt.figure(figsize=(8, 5))

    plt.scatter(moyenne, difference)

    plt.axhline(biais, linestyle="--", label=f"Biais = {biais:.2f}")
    plt.axhline(biais + 1.96 * sd, linestyle=":", label=f"± {1.96 * sd} SD")
    plt.axhline(biais - 1.96 * sd, linestyle=":")

    plt.xlabel("Moyenne des deux méthodes")
    plt.ylabel("Différence (code − NACCESS)")
    plt.title(f"Bland–Altman : comparaison des surfaces, {name} {molecule}")
    plt.legend()

    plt.show()

    plt.figure(figsize=(8, 5))
    # plt.plot(resultats["custom"].keys(), resultats["custom"].values(), label = "Custom")
    # plt.plot(resultats["naccess"].keys(), resultats["naccess"].values(), label = "Naccess")
    plt.fill_between(
        list(resultats["custom"].keys()),
        list(resultats["custom"].values()),
        label="Surface Custom",
        alpha = 0.4
    )

    plt.fill_between(
            list(resultats["custom"].keys()),
            list(resultats["naccess"].values()),
            label="Surface Naccess",
        alpha = 0.4
        )
    plt.legend()
    plt.xticks([1] + list(range(25, max(resultats["naccess"].keys()) + 1, 25)))
    plt.title(f"Comparaison des surfaces par résidus, {name} {molecule}")
    plt.xlabel("Numéro de résidu")
    plt.ylabel("Surface")
    plt.show()

    spearman, p_spearman = spearmanr(
        list(resultats["naccess"].values()),
        list(resultats["custom"].values())
    )

    pearson, p_pearson = pearsonr(
        list(resultats["naccess"].values()),
        list(resultats["custom"].values())
    )

    print(f"Spearman : {spearman:.4f}")
    print(f"Spearman p-value: {p_spearman:.4f}")
    print(f"Pearson  : {pearson:.4f}")
    print(f"Pearson p-value: {p_pearson:.4f}")