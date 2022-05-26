import os

from setuptools import setup
from src.version import __version__


def getdes():
    des = ""
    with open(os.path.join(os.getcwd(), "README.md")) as fi:
        des = fi.read()
    return des


setup(
    name="fsplit",
    version=__version__,
    packages=["fsplit"],
    package_dir={"fsplit": "src"},
    author="Deng Yong",
    author_email="yodeng@tju.edu.cn",
    url="",
    license="BSD",
    install_requires=[],
    python_requires='>=2.7.10, <=3.11',
    long_description=getdes(),
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'fsplit = fsplit.main:main',
        ]
    }
)
