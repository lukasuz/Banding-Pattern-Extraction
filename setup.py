from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="banding_pattern_extraction",
    version="0.1.0",
    author="Lukas Uzolas",
    author_email="lukas@uzolas.com",
    description="Package that allows for banding pattern extraction of metaphase chromosome images",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="lukas.uzolas.com",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)