from setuptools import setup

setup(
    name='cftAttachment',
    version='0.1.0',
    description='CFT Attachment calculator',
    url='https://edsaac.github.io',
    author='Edwin Saavedra C.',
    author_email='esaavedrac@u.northwestern.edu',
    license='MIT License',
    packages=['cftAttachment'],
    install_requires=['numpy',
                      'matplotlib',
                      'pandas',
                      'natsort',
                      'pyvista'
                      ],

    classifiers=['None'],
)