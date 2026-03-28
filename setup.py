from setuptools import setup, find_packages

setup(
    name="brca_harmonized",
    version="0.1.0",
    description="Cross-study BRCA harmonized database: TCGA-BRCA + METABRIC",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.5",
        "numpy>=1.23",
        "openpyxl>=3.0",
        "pyliftover>=0.4",
        "combat>=0.3.3",
    ],
    extras_require={
        "dev": ["pytest>=7.0"],
    },
)
