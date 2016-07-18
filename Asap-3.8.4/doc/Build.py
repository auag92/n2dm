#!/usr/bin/env python
from ASE.Utilities.BuildCAMPWebPages import Builder
from os.path import join, basename, split
from distutils.dep_util import newer

corner = """\
<a href="HOME/index.html">
  <img src="HOME/Asap-tmp-logo.png" align="top" width="133"
       height="55" alt="[ASAP]" border="0"/>
</a>"""

top = [
       ('ASAP introduction', 'intro.html'),
       ('documentation', 'documentation.html'),
       ('manual', 'manual/manual.html'),
       ('reference', 'reference.html'),
       ('examples', 'examples/examples.html'),
#       ('tutorials', 'tut/tutorials.html'),
#       ('faq', 'faq.html')
       ]

side = [('download', 'download.html'),
        ('mailing lists', 'lists.html'),
#        ('development', 'dev/development.html'),
#        ('publications', 'publications.html'),
#        ('technology', 'technology.html'),
#        ('reference', 'http://www.camp.dtu.dk/~schiotz/comp/Asap/epydoc/index.html'),
        ('bugs!', 'bugs.html')]

figures = ["""\
           <a href="http://www.fysik.dtu.dk">
             <img src="http://www.fysik.dtu.dk/images/template/logo-sml.gif"
                  align="top" width="125"
                  alt="[CAMP]" border="0"/>
           </a>""",
           """\
           <a href="http://validator.w3.org/check/referer">
             <img src="http://www.w3.org/Icons/valid-xhtml10"
                  alt="Valid XHTML 1.0!" height="31" width="88"
                  border="0"/>
           </a>""",
           """\
           <a href="http://www.fysik.dtu.dk/CAMPOS/ASE2">
             <img src="HOME/ASE-powered1.gif" align="top"
                  alt="[ASE-powered]" border="0"/>
           </a>""",
           """\
           <a href="http://www.python.org/">
             <img
               src="http://www.fysik.dtu.dk/campos/ASE/PythonPoweredSmall.gif"
               align="top" width="55" height="22"
               alt="[Python Powered]" border="0"/>
           </a>"""]


class ASAPBuilder(Builder):
    def Hook(self, dir, files):
        update = False
        if dir == './manual':
            reST = join(dir, 'manual.html')
            for f in files:
                if f.endswith('.txt') and newer(join(dir, f), reST):
                    update = True
                    break
            del files[:]
            files.append('manual.txt')
            
        return update


builder = ASAPBuilder(corner, top, side, figures)
builder.Build()
