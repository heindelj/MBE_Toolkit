
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='pyMD',  
     version='0.1',
     author="Joseph Heindel",
     author_email="heindelj@uw.edu",
     description="A molecular dynamics packine in python",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/heindelj/pyMD",
     install_requires=['numpy', 'tidynamics'],
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3"
     ],
 )
