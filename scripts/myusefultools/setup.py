from setuptools import setup

setup(
    name='myusefultools',
    version='0.1.0',
    description='Handle stuff from python',
    url='https://edsaac.github.io',
    author='Edwin Saavedra C.',
    author_email='esaavedrac@u.northwestern.edu',
    license='MIT License',
    packages=['myusefultools'],
    install_requires=['numpy',
                      'matplotlib',
                      'pandas',
                      'natsort'
                      ],

    classifiers=['None'],
)