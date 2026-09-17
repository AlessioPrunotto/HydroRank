from importlib.metadata import version

import hydrarank


def test_version_comes_from_installed_package_metadata():
    assert hydrarank.__version__ == version("hydrarank")
