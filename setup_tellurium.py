"""Setup with Tellurim."""

import setup_base

TELLURIUM = "tellurium"

install_requires = list(setup_base.INSTALL_REQUIRES)
install_requires.append(TELLURIUM)

setup_base.doSetup(install_requires)
