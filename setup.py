from setuptools import setup, find_packages

setup(
    name="intronxtract",
    version="1.0.0",
    packages=find_packages(),
    author="Carmen Lafuente Sanz",
    author_email="lafucarmen@gmail.com",
    description="Extract intron and flank sequences from BAM files",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/clafuent/intronxtract", 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CeCILL FREE SOFTWARE LICENSE AGREEMENT",
    ],
    python_requires='>=3.6',
    install_requires=[
        "pysam",
        "biopython",
    ],
    scripts=['scripts/intronxtract'] 
)
