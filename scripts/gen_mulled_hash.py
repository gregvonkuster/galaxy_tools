#!/usr/bin/env python3
# pip install galaxy-tool-util

import sys

from galaxy.tool_util.deps.mulled.util import build_target, v2_image_name


targets = []

for arg in sys.argv[1:]:
    assert '=' in arg
    name, version = arg.split('=')
    version.strip('=')  # in case you write ==
    targets.append(build_target(name, version))

name = v2_image_name(targets)

print(name)
