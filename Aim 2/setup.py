import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="corbus-python",
    version="0.0.1",
    author="David Randall Stokes",
    author_email="dstokes@mide.com",
    description="High-level API for Mide's CorBus haptics system",
    long_description=long_description,
    long_description_content_type="text/markdown; charset=UTF-8; variant=GFM",
    url="https://github.com/MideTechnology/corbus-python",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='mide corbus haptics api',
    python_requires='~=2.7',
    package_dir={'': '.'},
    package_data={
        '': ['hardware/*.xml', 'server/icons/*', 'tests/*.mkv']
    },
    install_requires=[
        'ebmlite==1.0.1',
        'numpy==1.16.2',
        'pyserial==3.4'
    ]
)
