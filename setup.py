from setuptools import setup, find_packages

setup(
    name='PyNetMHCpan',
    version='0.1.5',
    packages=find_packages(),
    url='https://github.com/kevinkovalchik/PyNetMHCpan',
    python_requires='>=3.7',
    license='MIT',
    author='Kevin Kovalchik',
    include_package_data=True,
    author_email='',
    description='A simple tool for using NetMHCpan and NetMHCIIpan using multiple CPUs in a Python environment.',
    install_requires=["numpy", "pandas"]
)
