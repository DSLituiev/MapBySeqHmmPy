from setuptools import setup, find_packages

version = "0.0.1"

setup(name="MapBySeqHmm",
    version=version,
    maintainer="Dmytro S Lituiev",
    maintainer_email="d.lituiev@gmail.com",
    license="http://opensource.org/licenses/BSD-3-Clause",
    platforms=["any"],
    long_description="""  """,

    packages=find_packages(exclude=['examples', 'tests']),

    include_package_data=True, # see MANIFEST.in
    zip_safe=False,
    entry_points={
	    'console_scripts':
	    [
	    'fm = MapBySeqHmm.MapBySeqHmmPy:filterMask',
	    ]
    },
)

