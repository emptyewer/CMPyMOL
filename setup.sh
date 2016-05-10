#!/usr/bin/env bash
rm -rf dist/*
rm -rf build/*
/usr/bin/python setup.py py2app --includes sip
echo "Unzipping Template DMG..."
bunzip2 -c template_dmg/template.dmg.bz2 > dist/temp.dmg
echo "Mount template_dmg/temp.dmg file, copy the DEEPN.app file to the mounted DEEPN folder, unmount temp.dmg, compress temp.dmg (using Disk Utility) and delete temp .dmg and DEEPN.app file."