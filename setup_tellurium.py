"""Setup with Tellurim."""

import setup

TELLURIUM = "tellurium"

install_requires = list(setup.INSTALL_REQUIRES)
install_requires.append(TELLURIUM)

setup.doSetup(install_requires)
