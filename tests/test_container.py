"""
Tests holding the container image's definition to the repository.

Dockerfile builds the image and installs requirements.lock. A package pinned in
requirements.txt must carry the same version in the lock, every locked package
must be pinned exactly with its hashes, the image must install from the lock and
from nothing else, and the paths holding data, restricted inputs or derived
outputs must stay out of the image. The tests read the files, so they run
without a container runtime.
"""

import pathlib
import re
import shlex

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "Dockerfile"
LOCK = ROOT / "requirements.lock"
REQUIREMENTS = ROOT / "requirements.txt"
IGNORE = ROOT / ".dockerignore"

# The paths the image is built from. Anything else reaching a COPY is a change
# to what the image holds and has to be made here as well.
# The image carries its own definition, so these tests also run inside it.
COPIED = {"Dockerfile", ".dockerignore", "requirements.txt", "requirements.lock",
          "LICENSE", "LICENSES.md", "CITATION.cff", "README.md", "MANIFEST.in",
          "pyproject.toml", "pytest.ini", "code/", "docs/", "micecoff/", "tests/"}
EXCLUDED = ["data/", "quarantine/", "outputs/", ".venv/", ".git/"]


def name(text):
    return re.sub(r"[-_.]+", "-", text).lower()


def instructions():
    """Instruction and arguments, with continuations joined and comments dropped."""
    found, current = [], ""
    for raw in DOCKERFILE.read_text().splitlines():
        line = raw.strip()
        if not current and (not line or line.startswith("#")):
            continue
        if line.endswith("\\"):
            current += line[:-1] + " "
            continue
        current += line
        keyword, _, arguments = current.partition(" ")
        found.append((keyword.upper(), arguments.strip()))
        current = ""
    assert not current, "the Dockerfile ends inside a continued instruction"
    return found


def requirement_pins():
    pins = {}
    for line in REQUIREMENTS.read_text().splitlines():
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s;]+)", line.strip())
        if match:
            pins[name(match.group(1))] = match.group(2)
    return pins


def locked():
    """Each locked package's version, marker and hashes."""
    entries = {}
    for block in re.split(r"\n(?=[A-Za-z0-9])", LOCK.read_text()):
        if block.startswith("#") or not block.strip():
            continue
        first = block.splitlines()[0].rstrip(" \\")
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s;]+)(?:\s*;\s*(.+))?", first)
        assert match, f"not an exact pin: {first}"
        key = name(match.group(1))
        assert key not in entries, f"{key} is locked twice"
        entries[key] = {"version": match.group(2), "marker": match.group(3),
                        "hashes": re.findall(r"--hash=sha256:([0-9a-f]+)", block)}
    return entries


def test_the_base_image_is_the_recorded_python_pinned_by_digest():
    bases = [arguments for keyword, arguments in instructions() if keyword == "FROM"]
    assert len(bases) == 1
    match = re.fullmatch(r"python:(\d+\.\d+\.\d+)-[a-z]+@sha256:[0-9a-f]{64}", bases[0])
    assert match, f"the base image is not pinned by digest: {bases[0]}"
    recorded = re.search(r"^# Python (\d+\.\d+\.\d+)$", REQUIREMENTS.read_text(), re.M)
    assert recorded and match.group(1) == recorded.group(1)


def test_every_pin_in_requirements_is_locked_at_the_same_version():
    pins, lock = requirement_pins(), locked()
    assert pins
    for package, version in pins.items():
        assert package in lock, f"{package} is pinned and not locked"
        assert lock[package]["version"] == version, package
        assert lock[package]["marker"] is None, f"{package} is locked conditionally"


def test_every_locked_package_carries_its_hashes():
    lock = locked()
    assert len(lock) > len(requirement_pins())
    for package, entry in lock.items():
        assert entry["hashes"], f"{package} has no hash"
        assert all(len(digest) == 64 for digest in entry["hashes"]), package


def test_the_lock_states_the_python_it_was_compiled_for():
    recorded = re.search(r"^# Python (\d+\.\d+)\.\d+$", REQUIREMENTS.read_text(), re.M)
    assert f"for Python {recorded.group(1)}." in " ".join(LOCK.read_text().split())


def pip_installs():
    commands = []
    for keyword, arguments in instructions():
        if keyword != "RUN":
            continue
        for command in re.split(r"\s*&&\s*", arguments):
            words = shlex.split(command)
            if words[:2] == ["pip", "install"]:
                commands.append(words[2:])
    return commands


def test_the_image_installs_from_the_lock_and_nothing_else():
    installs = pip_installs()
    assert len(installs) == 2
    packages, package = installs
    assert sorted(packages) == sorted(["--require-hashes", "--only-binary", ":all:",
                                       "--no-deps", "-r", "requirements.lock"])
    assert sorted(package) == sorted(["--no-deps", "--no-build-isolation",
                                      "--no-index", "."])
    runs = " ".join(arguments for keyword, arguments in instructions()
                    if keyword == "RUN")
    for fetcher in ("apt-get", "apt ", "apk ", "curl ", "wget ", "conda ", "uv "):
        assert fetcher not in runs, fetcher


def test_the_build_backend_is_locked_at_a_version_that_builds_the_package():
    tomllib = pytest.importorskip("tomllib")
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    floor = re.fullmatch(r"setuptools>=(\d+)", project["build-system"]["requires"][0])
    assert floor, project["build-system"]["requires"]
    version = locked()["setuptools"]["version"]
    assert int(version.split(".")[0]) >= int(floor.group(1))
    for dependency in project["project"]["dependencies"]:
        assert name(dependency.split(">=")[0]) in locked()


def copied_sources():
    sources = []
    for keyword, arguments in instructions():
        assert keyword != "ADD", "ADD can fetch and unpack; the image uses COPY"
        if keyword == "COPY":
            words = shlex.split(arguments)
            assert not any(word.startswith("--") for word in words), words
            sources.extend(words[:-1])
    return sources


def test_only_the_named_paths_enter_the_image():
    sources = copied_sources()
    assert set(sources) == COPIED
    for source in sources:
        assert (ROOT / source).exists(), f"{source} does not exist"
        assert not any(source.startswith(path) for path in EXCLUDED), source


def test_data_restricted_inputs_and_outputs_stay_out_of_the_build_context():
    lines = {line.strip() for line in IGNORE.read_text().splitlines()
             if line.strip() and not line.startswith("#")}
    for path in EXCLUDED:
        assert path in lines, f"{path} is not excluded from the build context"


def test_the_image_runs_the_test_suite_by_default():
    commands = [arguments for keyword, arguments in instructions() if keyword == "CMD"]
    assert commands == ['["pytest"]']


def test_the_readme_gives_the_images_paths():
    workdir = [arguments for keyword, arguments in instructions()
               if keyword == "WORKDIR"]
    assert len(workdir) == 1
    section = (ROOT / "README.md").read_text().split("\n## Container image\n")[1]
    section = section.split("\n## ")[0]
    assert "docker build" in section
    # A path inside the image, in prose or after the colon of a mount, which
    # excludes the host side written as $PWD/data.
    inside = re.findall(r"(?<![\w$])((?:/[\w.-]+)+)/(data|outputs|quarantine)\b",
                        section)
    assert {mount for _, mount in inside} >= {"data", "outputs"}
    for prefix, mount in inside:
        assert prefix == workdir[0], f"{prefix}/{mount} is not under {workdir[0]}"
