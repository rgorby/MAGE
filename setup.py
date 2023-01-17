import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kaipy",
    version="0.0.2",
    author="Center for Geospace Storms",
    author_email="eric.winter@jhuapl.edu",
    description="Postprocessing code for models in the kaiju software system",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/aplkaiju/kaiju/src/master/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    tests_require=['nose'],
    test_suite='nose.collector',
)
