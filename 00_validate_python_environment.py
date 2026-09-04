from importlib.metadata import PackageNotFoundError, version
import sys


EXPECTED_PYTHON = (3, 12, 10)
EXPECTED_PACKAGES = {
    "anndata": "0.12.10",
    "numpy": "2.4.3",
    "pandas": "2.3.3",
    "scipy": "1.17.1",
    "statsmodels": "0.14.6",
}


if sys.version_info[:3] != EXPECTED_PYTHON:
    found = ".".join(map(str, sys.version_info[:3]))
    expected = ".".join(map(str, EXPECTED_PYTHON))
    raise SystemExit(f"Python version mismatch: found {found}, expected {expected}.")

mismatches = []
for package, expected in EXPECTED_PACKAGES.items():
    try:
        found = version(package)
    except PackageNotFoundError:
        mismatches.append(f"{package}=missing expected={expected}")
        continue
    if found != expected:
        mismatches.append(f"{package}={found} expected={expected}")

if mismatches:
    raise SystemExit("Python package version mismatch: " + "; ".join(mismatches))

print("Python environment validation passed.")
