
from distutils.core import setup
import snp3merge

sys.path.append('googlemaps')
import googlemaps


setup(name='googlemaps',
      version='1.0',
      author='John Kleint',
      author_email='py-googlemaps-general@lists.sourceforge.net',
      url='http://sourceforge.net/projects/py-googlemaps/',
      download_url='https://sourceforge.net/projects/py-googlemaps/files/',
      description='Easy geocoding, reverse geocoding, driving directions, and local search in Python via Google.',
      long_description=googlemaps.GoogleMaps.__doc__,
      package_dir={'': 'googlemaps'},
      py_modules=['googlemaps'],
      provides=['googlemaps'],
      keywords='google maps local search ajax api geocode geocoding directions navigation json',
      license='Lesser Affero General Public License v3',
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Developers',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2',
                   'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
                   'License :: OSI Approved :: GNU Affero General Public License v3',
                   'Topic :: Internet',
                   'Topic :: Internet :: WWW/HTTP',
                   'Topic :: Scientific/Engineering :: GIS',
                  ],
     )
