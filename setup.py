import setuptools
import os


with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()

with open(
    os.path.join("computer_aided_design_for_optical_instruments", "version.py")
) as f:
    txt = f.read()
    last_line = txt.splitlines()[-1]
    version_string = last_line.split()[-1]
    version = version_string.strip("\"'")

setuptools.setup(
    name="computer_aided_design_for_optical_instruments",
    version=version,
    description=(
        "A python package for computer-aided design (CAD) for optical instruments"
    ),
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/computer_aided_design_for_optical_instruments",
    author="Sebastian Achim Mueller",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=[
        "computer_aided_design_for_optical_instruments",
        "computer_aided_design_for_optical_instruments.segmented_mirror",
        "computer_aided_design_for_optical_instruments.light_field_camera",
    ],
    package_data={"computer_aided_design_for_optical_instruments": []},
    install_requires=[
        "optic_object_wavefronts",
        "merlict",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
)
