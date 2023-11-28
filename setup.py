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
    description=("Create complex sceneries of optical instruments in merlict"),
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/optical_instruments_for_merlict",
    author="Sebastian Achim Mueller",
    author_email="Sebastian Achim Mueller@mail",
    packages=[
        "optical_instruments_for_merlict",
        "optical_instruments_for_merlict.segmented_mirror",
        "optical_instruments_for_merlict.light_field_camera",
    ],
    package_data={"optical_instruments_for_merlict": []},
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
