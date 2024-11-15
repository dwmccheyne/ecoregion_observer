from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="ecoregion_observer",  # Replace with your desired package name
    version="1.0.0",  # Initial version
    author="David McCheyne",
    author_email="davidmccheyne@gmail.com",
    description="A Python package for fetching, searching, and visualizing ecoregion observations",
    long_description=long_description,
    long_description_content_type="text/markdown",  # README.md format
    url="https://github.com/dwmccheyne/ecoregion_observer",  # Your GitHub repo URL
    packages=find_packages(),  # Automatically find sub-packages
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "requests>=2.25.1",
        "geopandas>=0.10.2",
        "pandas>=1.3.0",
        "matplotlib>=3.4.3",
        "shapely>=1.8.0",
    ],
    entry_points={
        "console_scripts": [
            "ecoregion-search-and-map=ecoregion_observer.ecoregion_search_and_map:main",
        ],
    },
    include_package_data=True,  # Include non-code files specified in MANIFEST.in
)
