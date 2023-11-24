import setuptools
import os


with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()

with open(os.path.join("optical_instruments_for_merlict", "version.py")) as f:
    txt = f.read()
    last_line = txt.splitlines()[-1]
    version_string = last_line.split()[-1]
    version = version_string.strip("\"'")

setuptools.setup(
    name="optical_instruments_for_merlict",
    version=version,
    description=("This is optical_instruments_for_merlict."),
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/optical_instruments_for_merlict",
    author="Sebastian Achim Mueller",
    author_email="Sebastian Achim Mueller@mail",
    packages=[
        "optical_instruments_for_merlict",
    ],
    package_data={"optical_instruments_for_merlict": []},
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
)
