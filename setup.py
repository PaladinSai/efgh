from setuptools import setup, find_packages

setup(
    name="efgh",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    author="Ming",
    description="A easy and fast GWAS tool in Python.",
    # install_requires=[],  # 之后可以补充依赖
)