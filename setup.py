from setuptools import setup

setup(
    name='PyNetMHCpan',
    version='0.1.2',
    packages=['PyNetMHCpan'],
    url='https://github.com/kevinkovalchik/PyNetMHCpan',
    license='MIT',
    author='Kevin Kovalchik',
    include_package_data=True,
    author_email='',
    description='A simple tool for using NetMHCpan and NetMHCIIpan using multiple CPUs in a Python environment.',
    install_requires=["numpy", "pandas"]
)
