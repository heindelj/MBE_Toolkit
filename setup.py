
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='py_MD',  
     version='0.3.2',
     author="Joseph Heindel",
     author_email="heindelj@uw.edu",
     description="A molecular dynamics package in python",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/heindelj/pyMD",
     download_url = 'https://github.com/heindelj/pyMD/archive/v_032.tar.gz',
     install_requires=['numpy', 'tidynamics', 'ase'],
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3"
     ],
 )
