# coding: utf-8

"""This module contains utility methods which are used to manage package versioning."""

from __future__ import unicode_literals
from datetime import datetime
from os.path import abspath, dirname
from subprocess import PIPE, Popen as Process


VERSION_MAPPING = {'alpha': 'a', 'beta': 'b', 'rc': 'c'}


def get_git_changeset():
    """Returns a numeric identifier of the latest Git changeset."""
    timestamp = Process('git log --pretty=format:%ct --quiet -1 HEAD', stdout=PIPE, stderr=PIPE, shell=True,
                        cwd=dirname(dirname(dirname(abspath(__file__)))), universal_newlines=True).communicate()[0]

    try:
        timestamp = datetime.utcfromtimestamp(int(timestamp))
    except ValueError:
        return str()

    return timestamp.strftime('%Y%m%d%H%M%S')


def get_version(version):
    """Returns a PEP 386-compliant version number obtained from the specified 'version' parameter."""
    main = ('{0}.{1}.{2}' if version[2] else '{0}.{1}').format(*version)

    if version[3] == 'alpha' and version[4] == 0:
        return "{0}.dev{1}".format(main, get_git_changeset())

    if version[3] != 'final':
        return "{0}{1}{2}".format(main, VERSION_MAPPING.get(version[3]), version[4])

    return main

