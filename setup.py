import setuptools


setuptools.setup(name='AISS',
version='0.1',
description='AI for Solid State Physics Applications',
url='#',
author='David Villarreal',
install_requires=['nequip','ase'],
author_email='david.leocadio94@gmail.com',
packages=setuptools.find_packages(include=["AISS","AISS.*"]),
zip_safe=False)
