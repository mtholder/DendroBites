__all__ = ['induced_matrix_and_tree']
from dendrobites.induced_matrix_and_tree import induced_matrix_and_tree

###############################################################################
## PACKAGE METADATA
import collections
from dendropy import homedir

version_info = collections.namedtuple("dendrobites_version_info",
                                      ["major",
                                       "minor",
                                       "micro",
                                       "releaselevel"])(major=0,
                                                        minor=0,
                                                        micro=0,
                                                        releaselevel="a")
__project__ = "DendroBites"
__version__ = ".".join(str(s) for s in version_info[:4] if s != "")
__author__ = "Jeet Sukumaran and Mark T. Holder"
__copyright__ = "Copyright 2015 Jeet Sukumaran and Mark T. Holder."
__citation__ = "Sukumaran, J and MT Holder. 2015. DendroBites python library"

def revision_description():
    __revision__ = _get_revision_object()
    if __revision__.is_available:
        revision_text = " ({})".format(__revision__)
    else:
        revision_text = ""
    return revision_text

def _get_revision_object():
    from dendropy.utility import vcsinfo
    __revision__ = vcsinfo.Revision(repo_path=homedir())
    return __revision__


def name():
    return "{} {}{}".format(__project__, __version__, revision_description())

def description(dest=None):
    import sys
    import site
    if dest is None:
        dest = sys.stdout
    fields = collections.OrderedDict()
    fields["DendroBites version"] = name()
    fields["DendroBites location"] = homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    try:
        fields["Python site packages"] = site.getsitepackages() #pylint: disable=E1103
    except:
        pass
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))
