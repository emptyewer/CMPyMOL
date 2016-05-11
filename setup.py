import os
import sys
import glob
from setuptools import setup
if sys.platform == 'darwin':
    import py2app
elif sys.platform == 'win32':
    import py2exe
    import matplotlib

def find_data_files(sources, targets, patterns):
    """Locates the specified data-files and returns the matches
    in a data_files compatible format.

    source is the root of the source data tree.
        Use '' or '.' for current directory.
    target is the root of the target data tree.
        Use '' or '.' for the distribution directory.
    patterns is a sequence of glob-patterns for the
        files you want to copy.
    """

    ret = {}
    for i, source in enumerate(sources):
        target = targets[i]
        if glob.has_magic(source) or glob.has_magic(target):
            raise ValueError("Magic not allowed in src, target")
        pattern = os.path.join(source, patterns[i])
        for filename in glob.glob(pattern):
            if os.path.isfile(filename):
                targetpath = os.path.join(target, os.path.relpath(filename, source))
                path = os.path.dirname(targetpath)
                ret.setdefault(path, []).append(filename)
    return sorted(ret.items())


APP = ['CMPyMOL_2.0.py']
INCLUDES = ['sip', 'PyQt4', 'time', 'sys', 'os', 'math', 'distutils.spawn', 're', 'subprocess', 'copy', 'xmlrpclib',
            'json', 'numbers', 'hashlib', 'decimal', 'thread', 'threading', 'matplotlib', 'itertools', 'operator',
            'matplotlib.backends.backend_qt4agg', 'matplotlib.colors', 'matplotlib.cm', 'matplotlib.figure',
            'matplotlib.patches', 'matplotlib.gridspec', 'numpy']
EXCLUDES = ['PyQt4.QtDesigner', 'PyQt4.QtNetwork', 'PyQt4.QtOpenGL', 'PyQt4.QtScript', 'PyQt4.QtSql', 'PyQt4.QtTest',
            'PyQt4.QtWebKit', 'PyQt4.QtXml', 'PyQt4.phonon', 'PIL']
OPTIONS = {'argv_emulation': True,
           'iconfile' : 'icon/Icon.icns',
           'plist': {'CFBundleGetInfoString': 'CMPyMOL',
                     'CFBundleIdentifier': 'edu.uiowa.krishnamani.cmpymol',
                     'CFBundleShortVersionString': '2.0',
                     'CFBundleName': 'CMPyMOL',
                     'CFBundleVersion': '2.0',
                     'NSHumanReadableCopyright': '(c) 2016 Venkatramanan Krishnamani'},
           'includes': INCLUDES,
           'excludes': EXCLUDES,
           }

if sys.platform == 'darwin':
    DATA_FILES = find_data_files(['ui', 'tools', 'overlays', 'graphs', 'functions'],
                                 ['ui', 'tools', 'overlays', 'graphs', 'functions'],
                                 ['*.ui', '*', '*.py', '*.py', '*.py'])
    setup(
        app=APP,
        name='CMPyMOL',
        options={'py2app': OPTIONS},
        setup_requires=['py2app'],
        author='Venkatramanan Krishnamani',
        data_files=DATA_FILES
    )
elif sys.platform == 'win32':
    DATA_FILES = find_data_files(['ui', 'tools', 'overlays', 'graphs', 'functions', 'dlls'],
                                 ['ui', 'tools', 'overlays', 'graphs', 'functions', ''],
                                 ['*.ui', '*', '*.py', '*.py', '*.py', '*'])
    origIsSystemDLL = py2exe.build_exe.isSystemDLL
    def isSystemDLL(pathname):
            if os.path.basename(pathname).lower() in ("msvcp71.dll", "dwmapi.dll", "'msvcp90.dll'"):
                    return 0
            return origIsSystemDLL(pathname)
    py2exe.build_exe.isSystemDLL = isSystemDLL
    setup(
        version='2.0',
        description='CMPyMOL',
        author='Venkatramanan Krishnamani',
        windows=[{"script": 'CMPyMOL_2.0.py',
                   "icon_resources": [(0, "icon/Icon.ico")],
                   "dest_base": "CMPyMOL"
                }],
        data_files=DATA_FILES + matplotlib.get_py2exe_datafiles(),
        options={"py2exe": {'includes': INCLUDES,
                            "optimize": 2,
                            "compressed": 2,
                            "bundle_files": 2,
                            "dist_dir": "dist\CMPyMOL",
                            "dll_excludes": ["numpy-atlas.dll"]
                            }},
        zipfile=None
    )
    setup(
        version='2.0',
        description='CMPyMOL',
        author='Venkatramanan Krishnamani',
        windows=[{"script": 'CMPyMOL_2.0.py',
                   "icon_resources": [(0, "icon/Icon.ico")],
                   "dest_base": "CMPyMOL"
                }],
        data_files=DATA_FILES + matplotlib.get_py2exe_datafiles(),
        options={"py2exe": {'includes': INCLUDES,
                            "optimize": 2,
                            "compressed": 2,
                            "bundle_files": 2,
                            "dist_dir": "dist\CMPyMOL",
                            "dll_excludes": ["numpy-atlas.dll"]
                            }},
        zipfile=None
    )
