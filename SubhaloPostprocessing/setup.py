from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="SubhaloPostprocessing",
    version="0.0.1",
    description="Classes for processing Subhalos",
    author="Erik Gillis",
    author_email="erik.gillis@mail.utoronto.ca",
    packages=find_packages(),
)
