from setuptools import setup

setup(
    name="muninn-biomass",
    version="1.0",
    description="Muninn extension for BIOMASS products",
    url="https://github.com/stcorp/muninn-biomass",
    author="S[&]T",
    license="BSD",
    py_modules=["muninn_biomass"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Environment :: Plugins",
    ],
    install_requires=["muninn"],
)
