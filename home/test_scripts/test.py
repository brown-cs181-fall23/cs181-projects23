import importlib

# DO NOT MODIFY THIS FILE #


# Create a function to check if a package is installed
def check_package(package_name):
    try:
        importlib.import_module(package_name)
        print(f"Package {package_name} is installed.")
        return True
    except ImportError:
        print(f"Package {package_name} is not installed.")
        return False


def test(packages_to_check):
    # Check each package in the list
    for package in packages_to_check:
        check_package(package)


if __name__ == "__main__":
    # List of packages to check
    packages_to_check = ["hmmlearn", "numpy", "pybedtools", "Bio", "pandas", "graphviz"]
    test(packages_to_check)
