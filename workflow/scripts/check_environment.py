"""Smoke-test the software required for the SAMap gut comparison."""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from shutil import which


PACKAGES = {
    "samap": "sc-samap",
    "scanpy": "scanpy",
    "anndata": "anndata",
    "numpy": "numpy",
    "pandas": "pandas",
    "snakemake": "snakemake",
}


def package_version(distribution: str) -> str:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "unknown"


def main() -> None:
    failures = []
    print("Python environment imports")
    for module_name, distribution in PACKAGES.items():
        try:
            import_module(module_name)
            print(f"  OK  {module_name}: {package_version(distribution)}")
        except Exception as exc:  # report binary/import incompatibilities too
            failures.append(f"{module_name}: {exc}")
            print(f"  FAIL {module_name}: {exc}")

    print("Command-line programs")
    for program in ("blastp", "makeblastdb"):
        path = which(program)
        if path:
            print(f"  OK  {program}: {path}")
        else:
            failures.append(f"{program}: not found on PATH")
            print(f"  FAIL {program}: not found on PATH")

    if failures:
        raise SystemExit("Environment check failed:\n- " + "\n- ".join(failures))

    print("Environment check passed.")


if __name__ == "__main__":
    main()
