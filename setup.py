try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from setuptools import find_packages

config = {
    'description': 'A set of tools to plan astronomical observations.',
    'author': 'iastro-pt',
    'url': 'https://github.com/iastro-pt/ObservationTools',
    'download_url': 'https://github.com/iastro-pt/ObservationTools',
    'author_email': 'daniel.andreasen@astro.up.pt',
    'version': '0.1',
    'license': 'MIT',
    'setup_requires': ['pytest-runner'],
    'tests_require': ['pytest', 'hypothesis'],
    # "PyAstronomy" when issue fixed.
    'install_requires': ["numpy", "astropy", "scipy", "matplotlib",
                         "astropy", "argparse", "ephem"],
    'extras_require': {
        'dev': ['check-manifest'],
        'tests': ['pytest', 'coverage', 'pytest-cov', 'python-coveralls', 'hypothesis'],
        'docs': ['sphinx >= 1.4'],
    },
    'packages': find_packages(exclude=['contrib', 'docs', 'tests', 'data']),
    'package_data': {
        # Inlcude the data files:
        '': ['data/*']},
    'scripts': ["visibility.py", "rv.py"],
    'name': 'ObservationTools',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    "classifiers": [
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English',
    ],
    # What does your project relate to?
    "keywords": ['Astronomy', 'Observation'],
}

setup(**config)
