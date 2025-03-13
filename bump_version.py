import re
import sys

def bump_version(file_path, pattern, new_version):
    with open(file_path, "r") as file:
        content = file.read()

    new_content = re.sub(pattern, new_version, content)

    with open(file_path, "w") as file:
        file.write(new_content)


def read_current_version(file_path, pattern):
    with open(file_path, "r") as file:
        content = file.read()

    match = re.search(pattern, content)
    if match:
        return match.group(0)
    else:
        raise ValueError(f"Version pattern not found in {file_path}")


def increment_version(current_version, part):
    major, minor, patch = map(int, re.findall(r"\d+", current_version))
    if part == "major":
        major += 1
        minor = 0
        patch = 0
    elif part == "minor":
        minor += 1
        patch = 0
    elif part == "patch":
        patch += 1
    else:
        raise ValueError("Invalid part to increment. Choose 'major', 'minor', or 'patch'.")
    return f"{major}.{minor}.{patch}"


def main(new_version_or_part):
    # Define file paths
    pyproject_file = "pyproject.toml"
    cmake_file = "CMakeLists.txt"
    pkgbuild_core_file = "packaging/PKGBUILD-core-avx2"
    pkgbuild_skyl_file = "packaging/PKGBUILD-skylake-avx512"
    makefilemanual_file = "src/makefile-manual"
    infoplist_file = "packaging/Info.plist"
    doxy_file = "docs/Doxyfile"

    # Define patterns to match the current version
    python_pattern = r'RASPA_VERSION\s*=\s*"\d+\.\d+\.\d+"'
    pyproject_pattern = r'version\s*=\s*"\d+\.\d+\.\d+"'
    cmake_pattern = r"VERSION \d+\.\d+\.\d+"
    pkgbuild_pattern = r"pkgver=\d+\.\d+\.\d+"
    makefilemanual_pattern = r"VERSION=\d+\.\d+\.\d+"
    infoplist_pattern = r"<string>\d+\.\d+\.\d+</string>"
    doxyfile_pattern = r"PROJECT_NUMBER\s*=\s*\d\.\d\.\d+"

    # Read current version from pyproject.toml
    current_version_match = re.search(r'"\d+\.\d+\.\d+"', read_current_version(pyproject_file, pyproject_pattern))
    if not current_version_match:
        raise ValueError("Current version not found in pyproject.toml")
    current_version = current_version_match.group(0).strip('"')

    # Determine the new version
    if new_version_or_part in ["major", "minor", "patch"]:
        new_version = increment_version(current_version, new_version_or_part)
    else:
        new_version = new_version_or_part

    # Define new version strings
    new_python_version = f'RASPA_VERSION = "{new_version}"'
    new_pyproject_version = f'version = "{new_version}"'
    new_cmake_version = f"VERSION {new_version}"
    new_pkgbuild_version = f"pkgver={new_version}"
    new_makefilemanual_version = f"VERSION={new_version}"
    new_infoplist_version = f"<string>{new_version}</string>"
    new_doxyfile_version = f"PROJECT_NUMBER         = {new_version}"
    

    # Bump versions in respective files
    bump_version(pyproject_file, pyproject_pattern, new_pyproject_version)
    bump_version(cmake_file, cmake_pattern, new_cmake_version)
    bump_version(pkgbuild_core_file, pkgbuild_pattern, new_pkgbuild_version)
    bump_version(pkgbuild_skyl_file, pkgbuild_pattern, new_pkgbuild_version)
    bump_version(makefilemanual_file, makefilemanual_pattern, new_makefilemanual_version)
    bump_version(infoplist_file, infoplist_pattern, new_infoplist_version)
    bump_version(doxy_file, doxyfile_pattern, new_doxyfile_version)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python bump_version.py <new_version_or_part>")
    elif (sys.argv[1] not in ["major", "minor", "patch"]) and (not bool(re.match(r'^\d+\.\d+\.\d+$', sys.argv[1]))):
        print("Usage: python bump_version.py <new_version_or_part>")
    else:
        main(sys.argv[1])
