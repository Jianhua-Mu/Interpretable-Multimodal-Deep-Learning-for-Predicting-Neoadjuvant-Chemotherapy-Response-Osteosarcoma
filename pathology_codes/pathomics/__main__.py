"""Show the available pathology-only command-line workflows."""


def main() -> None:
    print(
        "Pathomics commands:\n"
        "  python -m pathomics.generate_splits --help\n"
        "  python -m pathomics.train --help\n"
        "  python -m pathomics.evaluate --help\n"
        "  python -m pathomics.extract_patches --help\n"
        "  python -m pathomics.extract_features --help\n"
        "  python -m pathomics.generate_heatmaps --help"
    )


if __name__ == "__main__":
    main()
